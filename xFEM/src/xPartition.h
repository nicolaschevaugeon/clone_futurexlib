#ifndef XPARTITION_H
#define XPARTITION_H
#include <functional>
#include <set>

// AOMD
#include "mEntity.h"

// xfiles
#include "xEntityFilter.h"

namespace xfem
{
using xPartition = std::set<AOMD::mEntity*>;
using xGetPartition = std::function<void(AOMD::mEntity*, xPartition&, xEntityFilter)>;
}  // namespace xfem
#endif  // XPARTITION_H
