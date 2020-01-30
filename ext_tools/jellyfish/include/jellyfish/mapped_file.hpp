/*  This file is part of Jellyfish.

    This work is dual-licensed under 3-Clause BSD License or GPL 3.0.
    You can choose between one of them if you use this work.

`SPDX-License-Identifier: BSD-3-Clause OR  GPL-3.0`
*/

#ifndef __JELLYFISH_MAPPED_FILE_HPP__
#define __JELLYFISH_MAPPED_FILE_HPP__

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <vector>
#include <errno.h>
#include <iostream>

#include <jellyfish/err.hpp>

namespace jellyfish {
class mapped_file {
protected:
  std::string  _path;
  char        *_base, *_end;
  size_t       _length;

  void map_(int fd) {
    struct stat stat;
    if(fstat(fd, &stat) < 0)
      throw ErrorMMap(err::msg() << "Can't stat file '" << _path << "'" << err::no);

    _length = stat.st_size;
    _base = (char*)mmap(NULL, _length, PROT_READ, MAP_SHARED, fd, 0);
    if(_base == MAP_FAILED) {
      _base = 0;
      throw ErrorMMap(err::msg() << "Can't mmap file '" << _path << "'" << err::no);
    }
    _end = _base + _length;
  }

  void map_(const char *filename) {
    int fd = open(filename, O_RDONLY);
    if(fd < 0)
      throw ErrorMMap(err::msg() << "Can't open file '" << filename << "'" << err::no);
    map_(fd);
    close(fd);
  }


public:
  define_error_class(ErrorMMap);
  mapped_file() : _path(), _base(0), _end(0), _length(0) { }
  explicit mapped_file(const char *filename)
  : _path(filename), _base(0), _end(0), _length(0)
  {
    map_(filename);
  }
  explicit mapped_file(int fd)
  : _path(), _base(0), _end(0), _length()
  {
    map_(fd);
  }
  mapped_file(mapped_file&& rhs)
  : _path(std::move(rhs._path)), _base(rhs._base), _end(rhs._end),
    _length(rhs._length)
  {
    rhs._base = 0;
  }

  ~mapped_file() {
    unmap();
  }

  void map(const char* filename) {
    unmap();
    map_(filename);
  }

  void map(int fd) {
    unmap();
    map_(fd);
  }

  void unmap() {
    if(!_base)
      return;
    munmap(_base, _length);
    _path.clear();
    _base   = 0;
    _length = 0;
  }

  mapped_file& operator=(mapped_file&& rhs) {
    _path     = std::move(rhs._path);
    _base     = rhs._base;
    rhs._base = 0;
    _end      = rhs._end;
    _length   = rhs._length;
    return *this;
  }

  void swap(mapped_file& rhs) {
    std::swap(_path, rhs._path);
    std::swap(_base, rhs._base);
    std::swap(_end, rhs._end);
    std::swap(_length, rhs._length);
  }

  char *base() const { return _base; }
  char *end() const { return _end; }
  size_t length() const { return _length; }
  std::string path() const { return _path; }

  // No error checking here. Should I throw something?
  const mapped_file & will_need() const {
    madvise(_base, _length, MADV_WILLNEED);
    return *this;
  }
  const mapped_file & sequential() const {
    madvise(_base, _length, MADV_SEQUENTIAL);
    return *this;
  }
  const mapped_file & random() const {
    madvise(_base, _length, MADV_RANDOM);
    return *this;
  }
  const mapped_file & lock() const {
    if(mlock(_base, _length) < 0)
      throw ErrorMMap(err::msg() << "Can't lock map in memory" << err::no);
    return *this;
  }

  char load() const {
    const long    sz     = sysconf(_SC_PAGESIZE);
    // Do not optimize. Side effect is that every page is accessed and
    // should now be in cache.
    volatile char unused = 0;
    for(const char *w = _base; w < _base + _length; w += sz)
      unused ^= *w;
    return unused;
  }
};
inline void swap(mapped_file& a, mapped_file& b) { a.swap(b); }

// class mapped_files_t : public std::vector<mapped_file> {
// public:
//   mapped_files_t(int nb_files, char *argv[]) {
//     for(int j = 0; j < nb_files; j++)
//       push_back(mapped_file(argv[j]));
//   }

//   mapped_files_t(int nb_files, char *argv[], bool sequential) {
//     for(int j = 0; j < nb_files; j++) {
//       push_back(mapped_file(argv[j]));
//       if(sequential)
//         end()->sequential();
//     }
//   }
// };

// // File mapped on demand.
// class lazy_mapped_file_t : public mapped_file {
//   std::string       _path;
//   volatile bool     done;
//   volatile long     used_counter;

// public:
//   explicit lazy_mapped_file_t(const char *path) :
//     mapped_file((char *)0, (size_t)0),
//     _path(path), done(false), used_counter(0) {}

//   void map() {
//     used_counter = 1;
//     done = false;
//     mapped_file::map(_path.c_str());
//   }
//   void unmap() {
//     done = true;
//     dec();
//   }

//   void inc() {
//     atomic::gcc::fetch_add(&used_counter, (long)1);
//   }
//   void dec() {
//     long val = atomic::gcc::add_fetch(&used_counter, (long)-1);
//     if(done && val == 0)
//       mapped_file::unmap();
//   }
// };

// class lazy_mapped_files_t : public std::vector<lazy_mapped_file_t> {
// public:
//   lazy_mapped_files_t(int nb_files, char *argv[]) {
//     for(int j = 0; j < nb_files; j++)
//       push_back(lazy_mapped_file_t(argv[j]));
//   }

//   lazy_mapped_files_t(int nb_files, char *argv[], bool sequential) {
//     for(int j = 0; j < nb_files; j++) {
//       push_back(lazy_mapped_file_t(argv[j]));
//       if(sequential)
//         end()->sequential();
//     }
//   }
// };

}
#endif
