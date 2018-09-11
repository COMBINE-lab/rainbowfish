#include "rb-filesystem.hpp"
#include <sys/stat.h>


namespace rainbowfish {
    namespace fs {

        // Taken from http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
        bool FileExists(const char *path) {
            struct stat fileStat;
            if ( stat(path, &fileStat) ) {
                return false;
            }
            if ( !S_ISREG(fileStat.st_mode) ) {
                return false;
            }
            return true;
        }

        // Taken from http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
        bool DirExists(const char *path) {
            struct stat fileStat;
            if ( stat(path, &fileStat) ) {
                return false;
            }
            if ( !S_ISDIR(fileStat.st_mode) ) {
                return false;
            }
            return true;
        }

        void MakeDir(const char* path) {
            mkdir(path, ACCESSPERMS);
        }

    }
}
