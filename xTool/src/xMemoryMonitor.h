#ifndef _XMEMORYMONITOR_
#define _XMEMORYMONITOR_

/*
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 *          http://creativecommons.org/licenses/by/3.0/deed.en_US
 */

#if defined(_WIN32)
#include <psapi.h>
#include <windows.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <sys/resource.h>
#include <unistd.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif (defined(_AIX) || defined(__TOS__AIX__)) || \
    (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <cstdio>

#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif

/**
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 */
size_t getPeakRSS()
{
#if defined(_WIN32)
   /* Windows -------------------------------------------------- */
   PROCESS_MEMORY_COUNTERS info;
   GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
   return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) || \
    (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
   /* AIX and Solaris ------------------------------------------ */
   struct psinfo psinfo;
   int fd = -1;
   if ((fd = open("/proc/self/psinfo", O_RDONLY)) == -1) return (size_t)0L; /* Can't open? */
   if (read(fd, &psinfo, sizeof(psinfo)) != sizeof(psinfo))
   {
      close(fd);
      return (size_t)0L; /* Can't read? */
   }
   close(fd);
   return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
   /* BSD, Linux, and OSX -------------------------------------- */
   struct rusage rusage;
   getrusage(RUSAGE_SELF, &rusage);
#if defined(__APPLE__) && defined(__MACH__)
   return (size_t)rusage.ru_maxrss;
#else
   return (size_t)(rusage.ru_maxrss * 1024L);
#endif

#else
   /* Unknown OS ----------------------------------------------- */
   return (size_t)0L; /* Unsupported. */
#endif
}

/**
 * Returns the current resident set size (physical memory use) measured
 * in bytes, or zero if the value cannot be determined on this OS.
 */
size_t getCurrentRSS()
{
#if defined(_WIN32)
   /* Windows -------------------------------------------------- */
   PROCESS_MEMORY_COUNTERS info;
   GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
   return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
   /* OSX ------------------------------------------------------ */
   struct mach_task_basic_info info;
   mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
   if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &infoCount) != KERN_SUCCESS)
      return (size_t)0L; /* Can't access? */
   return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
   /* Linux ---------------------------------------------------- */
   long rss = 0L;
   FILE* fp = nullptr;
   if ((fp = fopen("/proc/self/statm", "r")) == nullptr) return (size_t)0L; /* Can't open? */
   if (fscanf(fp, "%*s%ld", &rss) != 1)
   {
      fclose(fp);
      return (size_t)0L; /* Can't read? */
   }
   fclose(fp);
   return (size_t)rss * (size_t)sysconf(_SC_PAGESIZE);

#else
   /* AIX, BSD, Solaris, and Unknown OS ------------------------ */
   return (size_t)0L; /* Unsupported. */
#endif
}

/* Usage

Simply call either function to get the current and peak resident set size, in bytes, of the current process:

size_t currentSize = getCurrentRSS( );
size_t peakSize    = getPeakRSS( );
*/

int parseLine(char* line)
{
   int i = strlen(line);
   while (*line < '0' || *line > '9') line++;
   line[i - 3] = '\0';
   i = atoi(line);
   return i;
}

/// Transform memory usage in Byte to memory in Human readable format (GB, KB etc ...)
std::string memoryHumanReadable(long int Byte)
{
   std::stringstream val;
   if (Byte / 1024 / 1024 / 1024)
   {
      val << Byte / 1024. / 1024. / 1024. << " GB";
      return val.str();
   }
   if (Byte / 1024 / 1024)
   {
      val << Byte / 1024. / 1024. << " MB";
      return val.str();
   }
   val << Byte << "  B";
   return val.str();
}

// return  system memory usage in Bytes Note That this is system dependent
int getSystemMemoryUsage()
{
   FILE* file = fopen("/proc/self/status", "r");
   int result = -1;
   char line[128];

   while (fgets(line, 128, file) != nullptr)
   {
      if (strncmp(line, "VmRSS:", 6) == 0)
      {
         result = parseLine(line);
         break;
      }
   }
   fclose(file);
   return result;
}

class xMemoryMonitor
{
  public:
   int start(std::string stage)
   {
      int index = strtoind[stage] = indtomemo.size();
      //    indtomemo.push_back(getSystemMemoryUsage());
      indtomemo.push_back(getCurrentRSS());
      return index;
   }
   void end(int id)
   {
      //    indtomemo[id] =getSystemMemoryUsage() -  indtomemo[id];
      indtomemo[id] = getCurrentRSS() - indtomemo[id];
   }

   void print(std::ostream& out = std::cout)
   {
      for (auto it = strtoind.begin(); it != strtoind.end(); ++it)
      {
         out.width(40);
         out.fill('_');
         out << std::left << it->first;
         out.fill(' ');
         out << std::right;
         out.width(10);
         out << memoryHumanReadable(indtomemo[it->second]) << std::endl;
      }
      out.width(40);
      out.fill('_');
      out << std::left;
      out << "Peak Memory ";
      out.width(10);
      out << memoryHumanReadable(getPeakRSS()) << std::endl;
   }

  private:
   std::map<std::string, int> strtoind;
   std::vector<int> indtomemo;
};
#endif
