import torch
import torch.nn as nn
import numpy as np
from scipy.stats import qmc
import matplotlib.pyplot as plt
import time
import os

# 设备配置
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# 可视化设置
try:
    plt.style.use('seaborn-v0_8')  # 使用新的seaborn样式名称
except:
    plt.style.use('ggplot')  # 如果seaborn-v0_8不可用，使用ggplot作为备选
plt.rcParams.update({
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'axes.unicode_minus': False,
    'figure.dpi': 150,
    'savefig.dpi': 300
})


class PhysicsInformedNN:
    def __init__(self, lambda_true=1.0, mu_true=0.5, Q=4.0,
                 lambda_init=0.8, mu_init=0.3,
                 hidden_layers=6, neurons=100,
                 n_samples=10000, use_alternating_optimization=False):
        # 材料参数
        self.lambda_true = lambda_true
        self.mu_true = mu_true
        self.Q = Q
        self.n_samples = n_samples
        self.use_alternating_optimization = use_alternating_optimization

        # 可训练参数
        self.lambda_param = nn.Parameter(torch.tensor(lambda_init, dtype=torch.float32, device=device))
        self.mu_param = nn.Parameter(torch.tensor(mu_init, dtype=torch.float32, device=device))

        # 神经网络
        self.displacement_net = self.build_network(input_dim=2, output_dim=2, hidden_layers=hidden_layers,
                                                   neurons=neurons)
        self.stress_net = self.build_network(input_dim=2, output_dim=3, hidden_layers=hidden_layers, neurons=neurons)

        # 训练数据（此时n_samples已定义）
        self.X, self.Y = self.create_lhs_samples()
        self.xy_train = torch.tensor(
            np.column_stack([self.X, self.Y]),
            dtype=torch.float32,
            device=device
        )

        # 预计算体力项（基于真实参数）
        self.f_x_true, self.f_y_true = self.precompute_body_forces()

        # 训练记录
        self.loss_history = []
        self.lambda_history = []
        self.mu_history = []

        # 优化器设置
        if self.use_alternating_optimization:
            # 为交替优化创建单独的优化器
            self.optimizer_lambda = torch.optim.Adam(
                list(self.displacement_net.parameters()) +
                list(self.stress_net.parameters()) +
                [self.lambda_param],
                lr=0.001
            )
            self.optimizer_mu = torch.optim.Adam(
                list(self.displacement_net.parameters()) +
                list(self.stress_net.parameters()) +
                [self.mu_param],
                lr=0.001
            )
        else:
            # 标准联合优化器
            self.optimizer = torch.optim.Adam(
                list(self.displacement_net.parameters()) +
                list(self.stress_net.parameters()) +
                [self.lambda_param, self.mu_param],
                lr=0.001
            )

    def precompute_body_forces(self):
        """基于真实材料参数预计算体力项"""
        pi = np.pi
        Q = self.Q

        # 使用真实参数计算体力
        f_x = self.lambda_true * (4 * pi ** 2 * np.cos(2 * pi * self.X) * np.sin(pi * self.Y) -
                                  pi * np.cos(pi * self.X) * Q * self.Y ** 3) + \
              self.mu_true * (9 * pi ** 2 * np.cos(2 * pi * self.X) * np.sin(pi * self.Y) -
                              pi * np.cos(pi * self.X) * Q * self.Y ** 3)

        f_y = self.lambda_true * (-3 * np.sin(pi * self.X) * Q * self.Y ** 2 +
                                  2 * pi ** 2 * np.sin(2 * pi * self.X) * np.cos(pi * self.Y)) + \
              self.mu_true * (-6 * np.sin(pi * self.X) * Q * self.Y ** 2 +
                              2 * pi ** 2 * np.sin(2 * pi * self.X) * np.cos(pi * self.Y) +
                              (pi ** 2 / 4) * np.sin(pi * self.X) * Q * self.Y ** 4)

        # 转换为PyTorch张量
        f_x_tensor = torch.tensor(f_x, dtype=torch.float32, device=device)
        f_y_tensor = torch.tensor(f_y, dtype=torch.float32, device=device)

        return f_x_tensor, f_y_tensor

    def create_lhs_samples(self):
        # 主体采样（强制转换为float32）
        main_samples = qmc.LatinHypercube(d=2).random(n=int(self.n_samples * 0.9)).astype(np.float32)
        main_samples = qmc.scale(main_samples, [0, 0], [1, 1])

        # 边界采样（强制转换为float32）
        n_edge = int(self.n_samples * 0.025)
        edges = np.concatenate([
            np.column_stack([np.random.rand(n_edge).astype(np.float32), np.zeros(n_edge, dtype=np.float32)]),
            np.column_stack([np.random.rand(n_edge).astype(np.float32), np.ones(n_edge, dtype=np.float32)]),
            np.column_stack([np.zeros(n_edge, dtype=np.float32), np.random.rand(n_edge).astype(np.float32)]),
            np.column_stack([np.ones(n_edge, dtype=np.float32), np.random.rand(n_edge).astype(np.float32)])
        ])

        return np.vstack([main_samples, edges])[:, 0], np.vstack([main_samples, edges])[:, 1]

    def build_network(self, input_dim, output_dim, hidden_layers, neurons):
        layers = [nn.Linear(input_dim, neurons), nn.Tanh()]  # 使用Tanh
        for _ in range(hidden_layers - 1):
            layers += [nn.Linear(neurons, neurons), nn.Tanh()]
        layers.append(nn.Linear(neurons, output_dim))
        return nn.Sequential(*layers).to(device)

    def exact_solution(self, x, y):
        """基于真实参数的解析解"""
        pi = np.pi
        Q = self.Q

        # 位移场
        u_x = np.cos(2 * pi * x) * np.sin(pi * y)
        u_y = 0.25 * Q * np.sin(pi * x) * y ** 4

        # 应变场
        eps_xx = -2 * pi * np.sin(2 * pi * x) * np.sin(pi * y)
        eps_yy = pi * Q * np.sin(pi * x) * y ** 3
        eps_xy = 0.5 * (pi * np.cos(2 * pi * x) * np.cos(pi * y) + Q * np.cos(pi * x) * y ** 3)

        # 应力场（使用真实材料参数）
        sigma_xx = self.lambda_true * (eps_xx + eps_yy) + 2 * self.mu_true * eps_xx
        sigma_yy = self.lambda_true * (eps_xx + eps_yy) + 2 * self.mu_true * eps_yy
        sigma_xy = 2 * self.mu_true * eps_xy

        return u_x, u_y, sigma_xx, sigma_yy, sigma_xy

    def compute_residuals(self, x, y):
        """计算物理残差（使用预计算的体力项）"""
        # 启用自动微分
        x.requires_grad_(True)
        y.requires_grad_(True)

        # 神经网络预测
        uv_pred = self.displacement_net(torch.stack([x, y], dim=1))
        u_pred = uv_pred[:, 0]  # x方向位移
        v_pred = uv_pred[:, 1]  # y方向位移

        sigma_pred = self.stress_net(torch.stack([x, y], dim=1))
        sigma_xx = sigma_pred[:, 0]
        sigma_yy = sigma_pred[:, 1]
        sigma_xy = sigma_pred[:, 2]

        # 计算应变分量 - 现在对每个分量单独计算梯度
        # x方向位移u对x的导数
        du_dx = torch.autograd.grad(u_pred.sum(), x, create_graph=True)[0]
        # y方向位移v对y的导数
        dv_dy = torch.autograd.grad(v_pred.sum(), y, create_graph=True)[0]
        # x方向位移u对y的导数
        du_dy = torch.autograd.grad(u_pred.sum(), y, create_graph=True)[0]
        # y方向位移v对x的导数
        dv_dx = torch.autograd.grad(v_pred.sum(), x, create_graph=True)[0]

        eps_xx = du_dx
        eps_yy = dv_dy
        eps_xy = 0.5 * (du_dy + dv_dx)

        # 计算应力导数 - 同样需要对每个分量单独计算
        dsigma_xx_dx = torch.autograd.grad(sigma_xx.sum(), x, create_graph=True)[0]
        dsigma_xy_dy = torch.autograd.grad(sigma_xy.sum(), y, create_graph=True)[0]
        dsigma_yy_dy = torch.autograd.grad(sigma_yy.sum(), y, create_graph=True)[0]
        dsigma_xy_dx = torch.autograd.grad(sigma_xy.sum(), x, create_graph=True)[0]

        # 使用预计算的体力项（固定值，不依赖待反演参数）
        f_x = self.f_x_true
        f_y = self.f_y_true

        # 平衡方程残差
        res_balance_x = dsigma_xx_dx + dsigma_xy_dy + f_x
        res_balance_y = dsigma_xy_dx + dsigma_yy_dy + f_y

        # 本构关系残差
        sigma_xx_constitutive = self.lambda_param * (eps_xx + eps_yy) + 2 * self.mu_param * eps_xx
        sigma_yy_constitutive = self.lambda_param * (eps_xx + eps_yy) + 2 * self.mu_param * eps_yy
        sigma_xy_constitutive = 2 * self.mu_param * eps_xy

        res_constitutive = torch.cat([
            (sigma_xx - sigma_xx_constitutive).unsqueeze(1),
            (sigma_yy - sigma_yy_constitutive).unsqueeze(1),
            (sigma_xy - sigma_xy_constitutive).unsqueeze(1)
        ], dim=1)

        return res_balance_x, res_balance_y, res_constitutive

    def loss_function(self):
        """计算损失函数（不执行优化步骤）"""
        # 获取训练数据坐标
        x = self.xy_train[:, 0]  # shape: (n_points,)
        y = self.xy_train[:, 1]  # shape: (n_points,)

        # --- 1. 数据匹配损失计算 ---
        # 计算解析解（注意转换为numpy数组后计算）
        with torch.no_grad():
            u_x_exact, u_y_exact, sigma_xx_exact, sigma_yy_exact, sigma_xy_exact = \
                self.exact_solution(x.detach().cpu().numpy(), y.detach().cpu().numpy())

            # 转换为PyTorch张量
            u_exact = torch.tensor(
                np.stack([u_x_exact, u_y_exact], axis=1),
                dtype=torch.float32, device=device
            )  # shape: (n_points, 2)

            sigma_exact = torch.tensor(
                np.stack([sigma_xx_exact, sigma_yy_exact, sigma_xy_exact], axis=1),
                dtype=torch.float32, device=device
            )  # shape: (n_points, 3)

        # 神经网络预测
        uv_pred = self.displacement_net(self.xy_train)  # shape: (n_points, 2)
        sigma_pred = self.stress_net(self.xy_train)  # shape: (n_points, 3)

        # 位移和应力的MSE损失
        mse_u = torch.mean(torch.sum((uv_pred - u_exact) ** 2, dim=1))  # 标量
        mse_sigma = torch.mean(torch.sum((sigma_pred - sigma_exact) ** 2, dim=1))  # 标量

        # --- 2. 物理约束损失计算 ---
        res_balance_x, res_balance_y, res_constitutive = self.compute_residuals(x, y)

        # 平衡方程残差 (动量守恒)
        mse_balance = torch.mean(res_balance_x ** 2 + res_balance_y ** 2)  # 标量

        # 本构关系残差 (应力-应变关系)
        mse_constitutive = torch.mean(torch.sum(res_constitutive ** 2, dim=1))  # 标量

        # --- 3. 参数约束（防止负值）---
        # 惩罚λ和μ为负的情况（relu(-x)在x<0时返回|x|）
        param_loss = torch.relu(-self.lambda_param) + torch.relu(-self.mu_param)  # 标量

        # --- 4. 组合损失函数 ---
        # 权重配置（需根据问题调整）
        # 动态调整权重（示例：物理约束权重随训练轮数增加）
        current_epoch = len(self.loss_history)
        balance_weight = min(1.0 + current_epoch / 5000, 5.0)  # 从1.0逐渐增至5.0
        constitutive_weight = min(1.0 + current_epoch / 5000, 5.0)

        loss_weights = {
            'data_u': 1.0,  # 位移数据匹配
            'data_sigma': 1.0,  # 应力数据匹配
            'balance': balance_weight,  # 平衡方程
            'constitutive': constitutive_weight,  # 本构关系
            'param_constraint': 0.1  # 参数约束
        }

        total_loss = (
                loss_weights['data_u'] * mse_u +
                loss_weights['data_sigma'] * mse_sigma +
                loss_weights['balance'] * mse_balance +
                loss_weights['constitutive'] * mse_constitutive +
                loss_weights['param_constraint'] * param_loss
        )

        return total_loss

    def train(self, n_epochs=10000):
        print("Starting training...")
        if self.use_alternating_optimization:
            print("Using alternating optimization strategy")
        else:
            print("Using joint optimization strategy")

        start_time = time.time()

        if self.use_alternating_optimization:
            # 为交替优化设置学习率调度器
            scheduler_lambda = torch.optim.lr_scheduler.ReduceLROnPlateau(
                self.optimizer_lambda, mode='min', factor=0.5, patience=200, verbose=False)
            scheduler_mu = torch.optim.lr_scheduler.ReduceLROnPlateau(
                self.optimizer_mu, mode='min', factor=0.5, patience=200, verbose=False)
        else:
            scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
                self.optimizer, mode='min', factor=0.5, patience=200, verbose=True)

        for epoch in range(n_epochs):
            if self.use_alternating_optimization:
                # 交替优化策略
                total_loss = self._alternating_optimization_step()

                # 更新学习率
                scheduler_lambda.step(total_loss)
                scheduler_mu.step(total_loss)
            else:
                # 标准联合优化
                self.optimizer.zero_grad()
                loss = self.loss_function()
                loss.backward()

                # 添加梯度裁剪
                torch.nn.utils.clip_grad_norm_(
                    list(self.displacement_net.parameters()) +
                    list(self.stress_net.parameters()) +
                    [self.lambda_param, self.mu_param],
                    max_norm=1.0
                )

                self.optimizer.step()
                scheduler.step(loss)
                total_loss = loss

            # 记录训练历史
            self.loss_history.append(total_loss.item())
            self.lambda_history.append(self.lambda_param.item())
            self.mu_history.append(self.mu_param.item())

            # 每1000轮重置采样点（可选）
            if epoch % 1000 == 0:
                self.X, self.Y = self.create_lhs_samples()
                self.xy_train = torch.tensor(
                    np.column_stack([self.X, self.Y]),
                    dtype=torch.float32,
                    device=device)
                # 重新计算体力项
                self.f_x_true, self.f_y_true = self.precompute_body_forces()

            # 打印进度
            if (epoch + 1) % 100 == 0:
                print(
                    f"Epoch {epoch + 1}/{n_epochs} | "
                    f"Loss: {total_loss.item():.3e} | "
                    f"Lambda: {self.lambda_param.item():.3f} (True: {self.lambda_true:.1f}) | "
                    f"Mu: {self.mu_param.item():.3f} (True: {self.mu_true:.1f})"
                )

        print(f"Training completed in {time.time() - start_time:.1f} seconds")

        # 保存结果
        os.makedirs("results_3_2", exist_ok=True)
        np.savez("results_3_2/training_history.npz",
                 loss=self.loss_history,
                 lambda_vals=self.lambda_history,
                 mu_vals=self.mu_history)

    def _alternating_optimization_step(self):
        """执行一步交替优化"""
        # 步骤1: 优化lambda和神经网络权重，冻结mu
        self.mu_param.requires_grad_(False)
        self.lambda_param.requires_grad_(True)

        self.optimizer_lambda.zero_grad()
        loss1 = self.loss_function()
        loss1.backward()

        # 梯度裁剪
        torch.nn.utils.clip_grad_norm_(
            list(self.displacement_net.parameters()) +
            list(self.stress_net.parameters()) +
            [self.lambda_param],
            max_norm=1.0
        )

        self.optimizer_lambda.step()

        # 步骤2: 优化mu和神经网络权重，冻结lambda
        self.lambda_param.requires_grad_(False)
        self.mu_param.requires_grad_(True)

        self.optimizer_mu.zero_grad()
        loss2 = self.loss_function()
        loss2.backward()

        # 梯度裁剪
        torch.nn.utils.clip_grad_norm_(
            list(self.displacement_net.parameters()) +
            list(self.stress_net.parameters()) +
            [self.mu_param],
            max_norm=1.0
        )

        self.optimizer_mu.step()

        # 恢复所有参数的梯度计算
        self.lambda_param.requires_grad_(True)
        self.mu_param.requires_grad_(True)

        # 返回最后一次的损失值
        return loss2

    def visualize_results(self):
        """可视化所有结果"""
        # 创建绘图目录
        os.makedirs("results_3_2/plots", exist_ok=True)

        # 1. 训练历史 - 损失函数曲线
        plt.figure(figsize=(8, 5))
        plt.semilogy(self.loss_history, linewidth=2)
        plt.xlabel("Iteration", fontsize=12)
        plt.ylabel("Loss Value (log scale)", fontsize=12)
        plt.title("Training Loss Curve", fontsize=14)
        plt.grid(False)
        plt.tight_layout()
        plt.close()

        # 2. 参数收敛曲线
        plt.figure(figsize=(10, 5))

        # Lambda收敛曲线
        plt.subplot(121)
        plt.plot(self.lambda_history, 'b-', linewidth=2, label="Estimated λ")
        plt.axhline(self.lambda_true, color='b', linestyle='--', linewidth=2, label="True λ")
        plt.xlabel("Iteration", fontsize=12)
        plt.ylabel("λ Value", fontsize=12)
        plt.title("λ Convergence", fontsize=14)
        plt.legend(fontsize=10)
        plt.grid(False)

        # Mu收敛曲线
        plt.subplot(122)
        plt.plot(self.mu_history, 'r-', linewidth=2, label="Estimated μ")
        plt.axhline(self.mu_true, color='r', linestyle='--', linewidth=2, label="True μ")
        plt.xlabel("Iteration", fontsize=12)
        plt.ylabel("μ Value", fontsize=12)
        plt.title("μ Convergence", fontsize=14)
        plt.legend(fontsize=10)
        plt.grid(False)

        plt.tight_layout()
        plt.close()

        # 3. 组合训练历史图
        plt.figure(figsize=(12, 5))
        plt.subplot(121)
        plt.semilogy(self.loss_history)
        plt.xlabel("Iteration")
        plt.ylabel("Loss")
        plt.title("Training Loss")
        plt.grid(False)

        plt.subplot(122)
        plt.plot(self.lambda_history, label=r"$\lambda$")
        plt.plot(self.mu_history, label=r"$\mu$")
        plt.axhline(self.lambda_true, color='C0', linestyle='--', label=r"True $\lambda$")
        plt.axhline(self.mu_true, color='C1', linestyle='--', label=r"True $\mu$")
        plt.xlabel("Iteration")
        plt.ylabel("Parameter Value")
        plt.title("Parameter Evolution")
        plt.legend()
        plt.grid(False)
        plt.tight_layout()
        plt.savefig("results_3_2/training_history.png")
        plt.close()

        # 预测结果（不再reshape）
        with torch.no_grad():
            uv_pred = self.displacement_net(self.xy_train)
            sigma_pred = self.stress_net(self.xy_train)
            u_pred = uv_pred[:, 0].cpu().numpy()  # 保持一维
            v_pred = uv_pred[:, 1].cpu().numpy()
            sigma_xx_pred = sigma_pred[:, 0].cpu().numpy()
            sigma_yy_pred = sigma_pred[:, 1].cpu().numpy()
            sigma_xy_pred = sigma_pred[:, 2].cpu().numpy()
            # 解析解（直接计算）
            u_x_exact, u_y_exact, sigma_xx_exact, sigma_yy_exact, sigma_xy_exact = \
                self.exact_solution(self.X, self.Y)
        # 位移场对比（分两次绘制）
        self._plot_comparison(u_pred, u_x_exact, "Displacement X", ["u_x"])
        self._plot_comparison(v_pred, u_y_exact, "Displacement Y", ["u_y"])
        # 应力场对比
        self._plot_comparison(sigma_xx_pred, sigma_xx_exact, "Stress XX", [r"$\sigma_{xx}$"])
        self._plot_comparison(sigma_yy_pred, sigma_yy_exact, "Stress YY", [r"$\sigma_{yy}$"])
        self._plot_comparison(sigma_xy_pred, sigma_xy_exact, "Stress XY", [r"$\sigma_{xy}$"])

        # 物理残差
        self._plot_residuals()

    def _plot_comparison(self, pred, exact, title, label):
        plt.figure(figsize=(10, 4))
        # 仅绘制单组对比图
        plt.scatter(self.X, self.Y, c=pred, cmap="jet", s=5, label="Predicted")
        plt.scatter(self.X, self.Y, c=exact, cmap="jet", s=5, alpha=0.5, label="Exact")
        plt.colorbar()
        plt.title(f"{title} Comparison")
        plt.legend()
        plt.savefig(f"results_3_2/plots/{title.lower()}_comparison.png")
        plt.close()

    def _plot_residuals(self):
        x = self.xy_train[:, 0].clone().requires_grad_(True)
        y = self.xy_train[:, 1].clone().requires_grad_(True)

        with torch.enable_grad():
            res_balance_x, res_balance_y, res_constitutive = self.compute_residuals(x, y)

        # 转换为NumPy并移除reshape
        residuals = {
            "Balance X": res_balance_x.detach().cpu().numpy(),
            "Balance Y": res_balance_y.detach().cpu().numpy(),
            "Constitutive XX": res_constitutive[:, 0].detach().cpu().numpy(),
            "Constitutive YY": res_constitutive[:, 1].detach().cpu().numpy(),
            "Constitutive XY": res_constitutive[:, 2].detach().cpu().numpy()
        }

        plt.figure(figsize=(15, 8))
        for i, (name, res) in enumerate(residuals.items()):
            plt.subplot(2, 3, i + 1)
            # 散点图可视化
            plt.scatter(self.X, self.Y, c=res, cmap="RdBu", s=2, vmin=-1, vmax=1)
            plt.colorbar()
            plt.title(name)
            plt.axis('off')

        plt.tight_layout()
        plt.savefig("results_3_2/plots/physics_residuals_scatter.png")
        plt.close()


if __name__ == "__main__":
    print("=== 测试标准联合优化 ===")
    # 初始化模型（标准优化）
    pinn_standard = PhysicsInformedNN(
        lambda_true=1.0, mu_true=0.5,
        lambda_init=0.8, mu_init=0.3,
        use_alternating_optimization=False
    )

    # 训练模型
    pinn_standard.train(n_epochs=5000)

    # 打印最终参数
    print("\nStandard Optimization Final Parameters:")
    print(f"  Lambda: {pinn_standard.lambda_param.item():.6f} (True: {pinn_standard.lambda_true})")
    print(f"  Mu:     {pinn_standard.mu_param.item():.6f} (True: {pinn_standard.mu_true})")

    print("\n=== 测试交替优化 ===")
    # 初始化模型（交替优化）
    pinn_alternating = PhysicsInformedNN(
        lambda_true=1.0, mu_true=0.5,
        lambda_init=0.8, mu_init=0.3,
        use_alternating_optimization=True
    )

    # 训练模型
    pinn_alternating.train(n_epochs=5000)

    # 可视化结果
    pinn_alternating.visualize_results()

    # 打印最终参数
    print("\nAlternating Optimization Final Parameters:")
    print(f"  Lambda: {pinn_alternating.lambda_param.item():.6f} (True: {pinn_alternating.lambda_true})")
    print(f"  Mu:     {pinn_alternating.mu_param.item():.6f} (True: {pinn_alternating.mu_true})")

    # 比较两种方法的结果
    print("\n=== 结果对比 ===")
    lambda_error_std = abs(pinn_standard.lambda_param.item() - pinn_standard.lambda_true)
    mu_error_std = abs(pinn_standard.mu_param.item() - pinn_standard.mu_true)
    lambda_error_alt = abs(pinn_alternating.lambda_param.item() - pinn_alternating.lambda_true)
    mu_error_alt = abs(pinn_alternating.mu_param.item() - pinn_alternating.mu_true)

    print(f"Standard Optimization Errors: λ={lambda_error_std:.6f}, μ={mu_error_std:.6f}")
    print(f"Alternating Optimization Errors: λ={lambda_error_alt:.6f}, μ={mu_error_alt:.6f}")
