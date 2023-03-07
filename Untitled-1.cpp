/* 最小二乘法拟合直线：AX+BY+C=0 */
bool xjLeastSquares::xjFitLineByLeastSquares(
std::vector<double> &xjParameters, const std::vector<xjPoint> &xjData)
{
	try
	{
		xjParameters.clear();
		if (xjData.size() < 2)
			return false;
 
		double A = 0.0, B = 0.0, C = 0.0;
 
		int xjPcount = xjData.size();
		double meanX = 0.0, meanY = 0.0;
		double sumXX = 0.0, sumXY = 0.0, sumYY = 0.0;
		for (int i = 0; i < xjPcount; i++)
		{
			meanX += xjData[i].x;
			meanY += xjData[i].y;
 
			sumXX += xjData[i].x * xjData[i].x;
			sumXY += xjData[i].x * xjData[i].y;
			sumYY += xjData[i].y * xjData[i].y;
		}
		meanX /= xjPcount;
		meanY /= xjPcount;
 
		sumXX -= xjPcount * meanX*meanX;
		sumXY -= xjPcount * meanX*meanY;
		sumYY -= xjPcount * meanY*meanY;
		if (abs(sumXX) < 1e-15)
		{
			A = 1.0;
			B = 0.0;
		}
		else
		{
			double ev = (sumXX + sumYY + sqrt((sumXX - sumYY)*(sumXX - sumYY) + 4 * sumXY*sumXY)) / 2.0;
			A = -sumXY;
			B = ev - sumYY;
			double norm = sqrt(A*A + B * B);
			A /= norm;
			B /= norm;
		}
 
		xjParameters.push_back(A);
		xjParameters.push_back(B);
		C = -(A * meanX + B * meanY);
		xjParameters.push_back(C);
 
		return true;
	}
	catch (const std::exception&)
	{
		return false;
	}
}