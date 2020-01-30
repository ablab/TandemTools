/*  This file is part of Jellyfish.

    This work is dual-licensed under 3-Clause BSD License or GPL 3.0.
    You can choose between one of them if you use this work.

`SPDX-License-Identifier: BSD-3-Clause OR  GPL-3.0`
*/

#ifndef __JELLYFISH_TEXT_DUMPER_HPP__
#define __JELLYFISH_TEXT_DUMPER_HPP__

#include <jellyfish/sorted_dumper.hpp>

namespace jellyfish {
template<typename Key, typename Val>
class text_writer {
public:
  void write(std::ostream& out, const Key& key, const Val val) {
   out << key << " " << val << "\n";
  }
};

template<typename storage_t>
class text_dumper : public sorted_dumper<text_dumper<storage_t>, storage_t> {
  typedef sorted_dumper<text_dumper<storage_t>, storage_t> super;
  text_writer<typename super::key_type, uint64_t> writer;

public:
  static const char* format;

  text_dumper(int nb_threads, const char* file_prefix, file_header* header = 0) :
    super(nb_threads, file_prefix, header)
  { }

  virtual void _dump(storage_t* ary) {
    if(super::header_) {
      super::header_->update_from_ary(*ary);
      super::header_->format(format);
    }
    super::_dump(ary);
  }

  void write_key_value_pair(std::ostream& out, typename super::heap_item item) {
    writer.write(out, item->key_, item->val_);
  }
};
template<typename storage_t>
const char* jellyfish::text_dumper<storage_t>::format = "text/sorted";

template<typename Key, typename Val>
class text_reader {
  std::istream& is_;
  char* buffer_;
  Key key_;
  Val val_;
  const RectangularBinaryMatrix m_;
  const size_t                  size_mask_;

public:
  text_reader(std::istream& is,
              file_header* header) :
    is_(is),
    buffer_(new char[header->key_len() / 2 + 1]),
    key_(header->key_len() / 2),
    m_(header->matrix()),
    size_mask_(header->size() - 1)
  { }

  const Key& key() const { return key_; }
  const Val& val() const { return val_; }
  size_t pos() const { return m_.times(key()) & size_mask_; }

  bool next() {
    is_ >> key_ >> val_;
    return is_.good();
  }
};
}

#endif /* __JELLYFISH_TEXT_DUMPER_HPP__ */
