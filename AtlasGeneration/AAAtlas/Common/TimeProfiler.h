
#ifndef __TIME_PROFILER_H
#define __TIME_PROFILER_H

#include <iostream>
#include <vector>
#include <map>
#include <time.h>
#include <assert.h>
#include <QString>


//#define CALCULATE_TIMING

class TimeProfiler
{
	public:
		TimeProfiler();
		void startedProcess(const QString& processID);
		void finishedProcess(const QString & processID, bool printTime = false);
		long getTotalProcessLengthInCycles(const QString & processID);
		QString getTotalProcessLengthInSeconds(const QString & processID);
		int getTotalProcessLengthInSecondsInt(const QString & processID);
		QString profile();
	
		void WriteProgress(const QString & process, int currProgress, int maxProgress);
		void WriteIntermediateTimes(std::vector<QString> & processes, int N);

#ifdef CALCULATE_TIMING
		static TimeProfiler * m_pGlobalSuperProfiler;
#endif

	protected:
		struct Process
		{
			Process();
			long m_iCurrStartTime;
//			std::vector<long> m_iRunningTimes;
			long m_iTotalRunningTime;
			QString m_sProcessName;
		};

		std::map<QString, Process> m_mIDToProcessMap;
};

#endif

