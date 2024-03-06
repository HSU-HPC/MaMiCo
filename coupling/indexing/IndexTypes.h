#pragma once

#include "coupling/indexing/CellIndex.h"

using I00 = coupling::indexing::CellIndex<3>;
using I01 = coupling::indexing::CellIndex<3,coupling::indexing::IndexTrait::vector>;
using I02 = coupling::indexing::CellIndex<3,coupling::indexing::IndexTrait::local>;
using I03 = coupling::indexing::CellIndex<3,coupling::indexing::IndexTrait::vector,coupling::indexing::IndexTrait::local>;
using I04 = coupling::indexing::CellIndex<3,coupling::indexing::IndexTrait::md2macro>;
using I05 = coupling::indexing::CellIndex<3,coupling::indexing::IndexTrait::vector,coupling::indexing::IndexTrait::md2macro>;
using I06 = coupling::indexing::CellIndex<3,coupling::indexing::IndexTrait::local,coupling::indexing::IndexTrait::md2macro>;
using I07 = coupling::indexing::CellIndex<3,coupling::indexing::IndexTrait::vector,coupling::indexing::IndexTrait::local,coupling::indexing::IndexTrait::md2macro>;
using I08 = coupling::indexing::CellIndex<3,coupling::indexing::IndexTrait::noGhost>;
using I09 = coupling::indexing::CellIndex<3,coupling::indexing::IndexTrait::vector,coupling::indexing::IndexTrait::noGhost>;
using I10 = coupling::indexing::CellIndex<3,coupling::indexing::IndexTrait::local,coupling::indexing::IndexTrait::noGhost>;
using I11 = coupling::indexing::CellIndex<3,coupling::indexing::IndexTrait::vector,coupling::indexing::IndexTrait::local,coupling::indexing::IndexTrait::noGhost>;
using I12 = coupling::indexing::CellIndex<3,coupling::indexing::IndexTrait::md2macro,coupling::indexing::IndexTrait::noGhost>;
using I13 = coupling::indexing::CellIndex<3,coupling::indexing::IndexTrait::vector,coupling::indexing::IndexTrait::md2macro,coupling::indexing::IndexTrait::noGhost>;
using I14 = coupling::indexing::CellIndex<3,coupling::indexing::IndexTrait::local,coupling::indexing::IndexTrait::md2macro,coupling::indexing::IndexTrait::noGhost>;
using I15 = coupling::indexing::CellIndex<3,coupling::indexing::IndexTrait::vector,coupling::indexing::IndexTrait::local,coupling::indexing::IndexTrait::md2macro,coupling::indexing::IndexTrait::noGhost>;
