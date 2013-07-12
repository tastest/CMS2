#include "Analysis/user.h"
#include "Core/AnalysisManager.h"
#include "Analysis/user.h"
#include "Services/logger.h"
// -----------------------------------------------------------------------------
// BuildTable
// -----------------------------------------------------------------------------
void AnalysisManager::BuildTable()
{
  Add(new user);
}
