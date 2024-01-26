#include "InfoOut.h"
#include <ctime>
#include <cstdarg>

InfoOut InfoManage;
InfoOut::InfoOut()
{
	screen = stdout;
	file = NULL;
	time_t curtime;
	time(&curtime); //ȡ��ǰʱ���ʾ1970��1��1�յ���ǰʱ�̵�����
	tm tm1;
	localtime_s(&tm1, &curtime); //ת��Ϊ����ʱ��
	sprintf_s(logfilename, 256, "%04d%02d%02d_%02d%02d.log", tm1.tm_year + 1900, tm1.tm_mon + 1, tm1.tm_mday, tm1.tm_hour, tm1.tm_min);
	//fopen_s(&file, logfilename, "a+");
}
InfoOut::~InfoOut()
{
	if (file)
	{
		fclose(file);
	}
	if (screen)
	{
		fclose(screen);
	}
}
string InfoOut::tostring(InfoType infotype)
{
	switch (infotype)
	{
	case error:
		return  "Error";
	case warning:
		return "Warning";
	case debug:
		return "Debug";
	case Info:
		return "Info";
	default:
		return NULL;
		break;
	}
}
void InfoOut::writeLogMsgtoScreen(InfoType infotype, const char* _fun, int _line, const char* format, ...)
{
	FILE* fp = screen;
	string str = tostring(infotype);
	fprintf_s(fp, "%s:(%d %s)", str.c_str(), _line, _fun);
	va_list args;
	va_start(args, format);
	vfprintf(fp, format, args);
	va_end(args);
	fprintf_s(fp, "\n");
}
void InfoOut::writeLogMsg(InfoType infotype, const char* format, ...)
{
	FILE* fp = screen;
	va_list args;
	va_start(args, format);
	vfprintf(fp, format, args);
	va_end(args);
	fprintf_s(fp, "\n");
}



