void LineFitLeastSquares(float *data_x, float *data_y, int data_n)  
{  
    float A = 0.0;  
    float B = 0.0;  
    float C = 0.0;  
    float D = 0.0;  
    float E = 0.0;  
    float F = 0.0;  
  
    for (int i=0; i<data_n; i++)  
    {  
        A += data_x[i] * data_x[i];  
        B += data_x[i];  
        C += data_x[i] * data_y[i];  
        D += data_y[i];  
    }  
  
    // 计算斜率a和截距b  
    float a, b, temp = 0;  
    if( temp = (data_n*A - B*B) )// 判断分母不为0  
    {  
        a = (data_n*C - B*D) / temp;  
        b = (A*D - B*C) / temp;  
    }  
    else  
    {  
        a = 1;  
        b = 0;  
    }  
  
    // 计算相关系数r  
    float Xmean, Ymean;  
    Xmean = B / data_n;  
    Ymean = D / data_n;  
  
    float tempSumXX = 0.0, tempSumYY = 0.0;  
    for (int i=0; i<data_n; i++)  
    {  
        tempSumXX += (data_x[i] - Xmean) * (data_x[i] - Xmean);  
        tempSumYY += (data_y[i] - Ymean) * (data_y[i] - Ymean);  
        E += (data_x[i] - Xmean) * (data_y[i] - Ymean);  
    }  
    F = sqrt(tempSumXX) * sqrt(tempSumYY);  
  
    float r;  
    r = E / F;  
}  