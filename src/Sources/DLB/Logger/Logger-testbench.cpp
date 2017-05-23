#include "Logger.h"
#include <iostream>

int main()
{
 
	TLogger::SetLevel(TLogger::BASIC);


	TLogger::Log(TLogger::BASIC, OUT_FMT_FIRST_SEPARATOR);
	TLogger::Log(TLogger::BASIC, OUT_FMT_CODE_NAME, "Logger hello, this is basic");
	TLogger::Log(TLogger::FULL, OUT_FMT_CODE_NAME, "Logger hello, this is full");
	TLogger::Log(TLogger::BASIC, OUT_FMT_SEPARATOR);
	TLogger::Log(TLogger::BASIC, OUT_FMT_LAST_SEPARATOR);
	TLogger::Log(TLogger::BASIC, OUT_FMT_END_OF_SIMULATION);
}