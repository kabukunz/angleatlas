
#include "TimeProfiler.h"
#ifdef CALCULATE_TIMING
	TimeProfiler * TimeProfiler::m_pGlobalSuperProfiler = new TimeProfiler();;
#endif

#include "QStringList.h"
TimeProfiler::Process::Process()
{
	m_iCurrStartTime = -1;
	m_iTotalRunningTime = 0;
}

TimeProfiler::TimeProfiler()
{
}

void TimeProfiler::startedProcess(const QString & processID)
{
	if (m_mIDToProcessMap[processID].m_iCurrStartTime!=-1)
		finishedProcess(processID);			//if process is currently running - finish it first
	
	m_mIDToProcessMap[processID].m_iCurrStartTime = clock();
	
	if (m_mIDToProcessMap[processID].m_sProcessName.isEmpty())
		m_mIDToProcessMap[processID].m_sProcessName = processID;
}

void TimeProfiler::finishedProcess(const QString& processID, bool printTime)
{
	Process & process = m_mIDToProcessMap[processID];
	long clockStamp = clock();
	if (process.m_iCurrStartTime==-1)
		process.m_iCurrStartTime = clockStamp;
	long diff = clockStamp - process.m_iCurrStartTime;
	process.m_iTotalRunningTime += (diff-1);
	process.m_iCurrStartTime = -1;

	if(printTime)
	{
		printf(">>>>>>> Finished: %s - %s\n", process.m_sProcessName.toStdString().c_str(), getTotalProcessLengthInSeconds(process.m_sProcessName).toStdString().c_str() );
	}
}

long TimeProfiler::getTotalProcessLengthInCycles(const QString& processID)
{
	return m_mIDToProcessMap[processID].m_iTotalRunningTime;
}

QString TimeProfiler::getTotalProcessLengthInSeconds(const QString& processID)
{
	long remainder = m_mIDToProcessMap[processID].m_iTotalRunningTime % CLOCKS_PER_SEC;
	remainder = (remainder * 1000) / CLOCKS_PER_SEC;
	QString retVal= QString::number(m_mIDToProcessMap[processID].m_iTotalRunningTime / CLOCKS_PER_SEC) + "." + QString::number(remainder) + "s";
	return retVal;

}

int TimeProfiler::getTotalProcessLengthInSecondsInt(const QString& processID)
{
	return m_mIDToProcessMap[processID].m_iTotalRunningTime / CLOCKS_PER_SEC;
}

QString TimeProfiler::profile()
{
	if (m_mIDToProcessMap.size()==0)
		return "Empty Time Profile";

	QString profile = "Time Profile:\n";
	for (std::map<QString, Process>::iterator iter = m_mIDToProcessMap.begin(); 
		 iter!=m_mIDToProcessMap.end(); iter++)
	{
		QString procName = iter->first;
		procName.replace(" ", "_");
		profile += QString("    ") + procName;
		profile += QString(" \t ") + getTotalProcessLengthInSeconds(iter->first) + "\n";
	}
	return profile;
}

void TimeProfiler::WriteProgress(const QString & process, int currProgress, int maxProgress)
{
	finishedProcess(process);	
	int TICK_SLOTS = 40;
	int currTicks = TICK_SLOTS * currProgress / maxProgress;
	std::cout<<"\r[";
	for (int i=0; i<TICK_SLOTS; i++)
	{
		if (i < currTicks)
			std::cout<<"+"<<std::flush;
		else
			std::cout<<" "<<std::flush;
	}
	
	int currTime = getTotalProcessLengthInSecondsInt(process);
	int remaining = (maxProgress-currProgress);
	int remainingTime = ((currProgress==0) ? (-1) : (currTime * remaining / currProgress));
	std::cout<<"] N="<<remaining<<" T={"<<currTime<<"s | "<<remainingTime<<"s}    "<<std::flush;
	startedProcess(process);	
}

void TimeProfiler::WriteIntermediateTimes(std::vector<QString> & processes, int N)
{
	std::cout<<"\r[";	
	for (int i=0; i<processes.size(); i++)
	{
		finishedProcess(processes[i]);
		std::cout<<(getTotalProcessLengthInSecondsInt(processes[i]) / N)<<"s\t"<<std::flush;
		startedProcess(processes[i]);	
	}
		
	std::cout<<"] / "<<N<<std::flush;
}







