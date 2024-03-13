import joblib
import xgboost as xgb
import seaborn as sns
import numpy as np
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.model_selection import cross_val_score
import pandas as pd
# from sklearn.experimental import enable_hist_gradient_boosting  # noqa
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor, ExtraTreesRegressor, AdaBoostRegressor, \
     HistGradientBoostingRegressor, VotingRegressor
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.datasets import load_iris, load_digits, load_diabetes, load_boston, load_wine, load_breast_cancer
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn import svm
from sklearn import tree
from sklearn import metrics
from sklearn.linear_model import PassiveAggressiveClassifier
from sklearn.preprocessing import FunctionTransformer, StandardScaler, MinMaxScaler
from sklearn.decomposition import PCA
from sklearn import linear_model
from sklearn.preprocessing import PolynomialFeatures
# import tensorflow as tf
from sklearn import datasets, linear_model
from sklearn.model_selection import cross_val_score, cross_validate
from sklearn.metrics import r2_score
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from scipy import stats
import time
import ggplot.datasets as gds
from scipy.optimize import minimize


def de_mean(x):
    x_bar = np.mean(x)
    return [x_i - x_bar for x_i in x]


def dot(v, w):
    return sum(v_i * w_i for v_i, w_i in zip(v, w))


def sum_of_squares(v):
    return dot(v, v)


def covariance(x, y):
    n = len(x)
    return dot(de_mean(x), de_mean(y)) / (n-1)


def correlation(x, y):
    stdev_x = np.std(x)
    stdev_y = np.std(y)
    if stdev_x > 0 and stdev_y > 0:
        return covariance(x, y) / stdev_x / stdev_y
    else:
        return 0


def rps_ams_auto_modeling(x, y):
     rps_ams = joblib.load("RPS-AMS_RankingsPredictor.m")
     x = np.array(x)
     y = np.array(y)
     y = np.ravel(y)

     np.random.seed(42)
     np.random.shuffle(x)
     np.random.seed(42)
     np.random.shuffle(y)

     x = x.astype(float)
     y = y.astype(float)

     x_average = np.mean(x, axis=0)
     y_average = np.mean(y, axis=0)

     col0 = x.shape[1]
     row0 = x.shape[0]
     col = 50
     x_average_add = np.zeros([col - col0], dtype=float)
     x_average = np.concatenate([x_average, x_average_add])

     n_split = 1
     # the number of times of test
     index = 1

     meta_feature_set = np.empty(((int(n_split * index)), int(60)), dtype=float)
     batch_size = int(row0 / n_split)

     # save the R2_score of all surrogate models
     target = np.empty(int(12), dtype=float)
     # Calculate all R2_score

     # construct the meta-feature dataset
     start = time.time()
     for k in range(index):
          y_ = y
          y_ = np.ravel(y_)
          y_ = np.array(y_)

          x_train = x
          y_train = y_

          x_train = np.array_split(x_train, n_split)
          y_train = np.array_split(y_train, n_split)

          for i in range(int(n_split)):
               # define the representative surrogate model
               svr_rbf = svm.SVR()
               gbr = GradientBoostingRegressor()
               et = ExtraTreesRegressor()

               # define the Polynomial Regression to calculate the linearity
               poly_Feature = PolynomialFeatures()
               linear = linear_model.LinearRegression()
               row = x_train[i].shape[0]
               # if necessary, apply the zero-padding method
               x_train_add = np.zeros([row, col - col0], dtype=float)
               x_train_reshape = np.concatenate([x_train[i], x_train_add], axis=1)
               # calculate the statistics moments
               y_train_average = np.mean(y_train[i], axis=0)
               y_train_std = np.std(y_train[i], axis=0)
               y_train_skewness = stats.skew(y_train[i])
               y_train_kurtosis = stats.kurtosis(y_train[i])
               # calculate the correlation coefficient
               for j in range(50):
                    meta_feature_set[k * (int(n_split)) + i, j] = correlation(x_train_reshape[:, j], y_train[i])

               meta_feature_set[k * (int(n_split)) + i, 50] = y_train_average - y_average
               meta_feature_set[k * (int(n_split)) + i, 51] = y_train_std
               meta_feature_set[k * (int(n_split)) + i, 52] = y_train_skewness
               meta_feature_set[k * (int(n_split)) + i, 53] = y_train_kurtosis
               meta_feature_set[k * (int(n_split)) + i, 54] = batch_size
               meta_feature_set[k * (int(n_split)) + i, 55] = col0
               # if R2 is at the exceptional value, set it as 0 or -1
               error_value = -1
               cv = 10
               x_train_batch, x_test_batch, y_train_batch, y_test_batch = train_test_split(x_train[i], y_train[i],
                                                                                           test_size=0.1,
                                                                                           random_state=42)
               # calculate the accuracy of the representative surrogate models
               et.fit(x_train_batch, y_train_batch)
               score_ET = r2_score(y_test_batch, et.predict(x_test_batch))
               if score_ET <= error_value:
                    score_ET = error_value
               meta_feature_set[k * (int(n_split)) + i, 56] = score_ET

               gbr.fit(x_train_batch, y_train_batch)
               score_GBR = r2_score(y_test_batch, gbr.predict(x_test_batch))
               if score_GBR <= error_value:
                    score_GBR = error_value
               meta_feature_set[k * (int(n_split)) + i, 57] = score_GBR

               linear.fit(x_train_batch, y_train_batch)
               score_Poly = r2_score(y_test_batch, linear.predict(x_test_batch))
               if score_Poly <= error_value:
                    score_Poly = error_value

               meta_feature_set[k * (int(n_split)) + i, 58] = score_Poly

               svr_rbf.fit(x_train_batch, y_train_batch)
               score_SVR_rbf = r2_score(y_test_batch, svr_rbf.predict(x_test_batch))
               if score_SVR_rbf <= error_value:
                    score_SVR_rbf = error_value
               meta_feature_set[k * (int(n_split)) + i, 59] = score_SVR_rbf

     num_algorithm = 12
     x_test = meta_feature_set

     num_row_test = np.shape(x_test)[0]
     # predict the R2 score of all surrogate model
     y_pred_rf = rps_ams.predict(x_test)
     y_pred_rf = np.reshape(y_pred_rf, (num_row_test, num_algorithm))
     # calculate the rankings of all surrogate models
     sort_predict = np.argsort(np.argsort(-np.mean(y_pred_rf, axis=0)))

     pipe_dic = {0: 'Random Forest', 1: 'Decision Tree', 2: 'extra tree', 3: 'gradient boosting', 4: 'ada boost ',
                 5: 'Hist Gradient boosting', 6: 'SVR kernel rbf',
                 7: 'SVR kernel poly', 8: 'SVR kernel sigmoid', 9: 'NuSVR kernel rbf',
                 10: 'NuSVR kernel poly', 11: 'NuSVR kernel sigmoid'}

     first_model_id = np.argwhere(sort_predict == 0)
     second_model_id = np.argwhere(sort_predict == 1)
     third_model_id = np.argwhere(sort_predict == 2)
     forth_model_id = np.argwhere(sort_predict == 3)
     cv = 3

     # model selection

     def selecting_model(model_id):
          best_model = 0
          best_model_score = 0
          best_model_r2_score = 0
          best_model_r2_std = 0
          if model_id == 0:
               best_model = RandomForestRegressor()
               best_model_score = cross_validate(best_model, x, y, cv=cv, scoring='r2', return_estimator=True)
               best_model_r2_score = best_model_score['test_score'].mean()
               best_model_r2_std = best_model_score['test_score'].std()
          elif model_id == 1:
               best_model = tree.DecisionTreeRegressor()
               best_model_score = cross_validate(best_model, x, y, cv=cv, scoring='r2', return_estimator=True)
               best_model_r2_score = best_model_score['test_score'].mean()
               best_model_r2_std = best_model_score['test_score'].std()
          elif model_id == 2:
               best_model = ExtraTreesRegressor()
               best_model_score = cross_validate(best_model, x, y, cv=cv, scoring='r2', return_estimator=True)
               best_model_r2_score = best_model_score['test_score'].mean()
               best_model_r2_std = best_model_score['test_score'].std()
          elif model_id == 3:
               best_model = GradientBoostingRegressor()
               best_model_score = cross_validate(best_model, x, y, cv=cv, scoring='r2', return_estimator=True)
               best_model_r2_score = best_model_score['test_score'].mean()
               best_model_r2_std = best_model_score['test_score'].std()
          elif model_id == 4:
               best_model = AdaBoostRegressor()
               best_model_score = cross_validate(best_model, x, y, cv=cv, scoring='r2', return_estimator=True)
               best_model_r2_score = best_model_score['test_score'].mean()
               best_model_r2_std = best_model_score['test_score'].std()
          elif model_id == 5:
               best_model = HistGradientBoostingRegressor()
               best_model_score = cross_validate(best_model, x, y, cv=cv, scoring='r2', return_estimator=True)
               best_model_r2_score = best_model_score['test_score'].mean()
               best_model_r2_std = best_model_score['test_score'].std()
          elif model_id == 6:
               best_model = svm.SVR(kernel='rbf')
               best_model_score = cross_validate(best_model, x, y, cv=cv, scoring='r2', return_estimator=True)
               best_model_r2_score = best_model_score['test_score'].mean()
               best_model_r2_std = best_model_score['test_score'].std()
          elif model_id == 7:
               best_model = svm.SVR(kernel='poly')
               best_model_score = cross_validate(best_model, x, y, cv=cv, scoring='r2', return_estimator=True)
               best_model_r2_score = best_model_score['test_score'].mean()
               best_model_r2_std = best_model_score['test_score'].std()
          elif model_id == 8:
               best_model = svm.SVR(kernel='sigmoid')
               best_model_score = cross_validate(best_model, x, y, cv=cv, scoring='r2', return_estimator=True)
               best_model_r2_score = best_model_score['test_score'].mean()
               best_model_r2_std = best_model_score['test_score'].std()
          elif model_id == 9:
               best_model = svm.NuSVR(kernel='rbf')
               best_model_score = cross_validate(best_model, x, y, cv=cv, scoring='r2', return_estimator=True)
               best_model_r2_score = best_model_score['test_score'].mean()
               best_model_r2_std = best_model_score['test_score'].std()
          elif model_id == 10:
               best_model = svm.NuSVR(kernel='poly')
               best_model_score = cross_validate(best_model, x, y, cv=cv, scoring='r2', return_estimator=True)
               best_model_r2_score = best_model_score['test_score'].mean()
               best_model_r2_std = best_model_score['test_score'].std()
          elif model_id == 11:
               best_model = svm.NuSVR(kernel='sigmoid')
               best_model_score = cross_validate(best_model, x, y, cv=cv, scoring='r2', return_estimator=True)
               best_model_r2_score = best_model_score['test_score'].mean()
               best_model_r2_std = best_model_score['test_score'].std()
          return best_model, best_model_r2_score, best_model_r2_std

     # calculate the types, accuracy and std of surrogate
     first_model, first_model_r2_score, first_model_r2_std = selecting_model(first_model_id)
     second_model, second_model_r2_score, second_model_r2_std = selecting_model(second_model_id)
     third_model, third_model_r2_score, third_model_r2_std = selecting_model(third_model_id)
     forth_model, forth_model_r2_score, forth_model_r2_std = selecting_model(forth_model_id)

     score_candidate = [first_model_r2_score, second_model_r2_score, third_model_r2_score, forth_model_r2_score]
     std_candidate = [first_model_r2_std, second_model_r2_std, third_model_r2_std, forth_model_r2_std]

     score = np.array(score_candidate)
     std_candidate = np.array(std_candidate)
     best_value = score.max()

     threshold = 0.95 * best_value

     r2_square = np.where(score > threshold, score, 0)
     r2_std = std_candidate

     # figure the weights

     cons = ({'type': 'eq', 'fun': lambda x: x[0] + x[1] + x[2] + x[3] - 1})
     bounds = [(0, 1)] * 4
     prediction_weights = np.random.random(4)

     def optimize_function(weights):
          res = 5 * weights * r2_std - 1 * weights * r2_square
          res = res.sum()
          return res

     res = minimize(optimize_function, prediction_weights, method='SLSQP', bounds=bounds, constraints=cons)
     end = time.time()

     # Ensemble the surrogate models
     clf = VotingRegressor(estimators=[('first', first_model), ('second', second_model),
                                       ('third', third_model), ('forth', forth_model)], weights=res.x)
     cv2 = 10
     # modeling
     clf_score = cross_validate(clf, x, y, cv=cv2, scoring='r2')
     print('time elapse:', end - start)
     print('R2 and std :', clf_score['test_score'].mean(), clf_score['test_score'].std())
     # display the candidate surrogate
     print(first_model, second_model, third_model, forth_model)
     return clf_score['test_score'].mean()