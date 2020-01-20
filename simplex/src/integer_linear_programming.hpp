#pragma once

#include "linear_programming.hpp"

namespace simplex {

  // ILP in SEF. Not for mixed linear programming
  template <typename FieldType>
  using ILPProblem = LPProblem<FieldType>;

  template <typename FieldType>
  bool result_is_integral(const Result<FieldType>&);

  template <typename FieldType>
  Result<FieldType>
  ip_solve(ILPProblem<FieldType>);

  template <typename FieldType>
  Result<FieldType>
  ip_solve_branch_and_cut(ILPProblem<FieldType>);

}

