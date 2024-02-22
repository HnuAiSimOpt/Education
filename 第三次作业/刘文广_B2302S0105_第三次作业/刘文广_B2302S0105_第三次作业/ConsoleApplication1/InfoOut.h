#ifndef INFOOUT_H
#define INFOOUT_H
#include <cstdio>
#include <string>

using std::string;
class InfoOut
{
private:
	FILE *file, *screen;
	char logfilename[256];
public:
	//构造函数
	enum InfoType
	{
		error,
		warning,
		debug,
		Info,
	};
	InfoOut();
	~InfoOut();
	void writeLogMsgtoScreen(InfoType infotype, const char* _fun, int _line, const char* format, ...);
	void writeLogMsg(InfoType infotype, const char* fomat, ...);
	string tostring(InfoType infotype);

};
extern InfoOut InfoManage;
/* 
 * ERROR: 报告错误但会退出
 */
#define ERROR(...) InfoManage.writeLogMsgtoScreen(InfoOut::error, __FUNCTION__, __LINE__, __VA_ARGS__); exit(1);
#define WARNING(...) InfoManage.writeLogMsgtoScreen(InfoOut::warning, __FUNCTION__, __LINE__, __VA_ARGS__)
#define DEBUG(...) InfoManage.writeLogMsgtoScreen(InfoOut::debug, __FUNCTION__, __LINE__, __VA_ARGS__)
#define INFO(...) InfoManage.writeLogMsg(InfoOut::Info, __VA_ARGS__)

#endif