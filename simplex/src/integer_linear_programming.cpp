#include "integer_linear_programming.hpp"
#include <queue>

namespace simplex {
  template <typename FieldType>
  Result<FieldType>
  solve(ILPProblem<FieldType> ilp)
  { return solve_branch_and_cut(std::move(ilp)); }

  template <typename FieldType>
  Result<FieldType>
  solve_branch_and_cut(ILPProblem<FieldType> ilp)
  {
    std::queue<ILPProblem<FieldType>> q;
    q.push(std::move(ilp));
    Result<FieldType> res;
    FieldType best = Eigen::NumTraits<FieldType>::lowest();
    while (!q.empty()) {
      ILPProblem<FieldType> top = q.front();
      q.pop();
      do {
        Result<FieldType> relax_res = lp_solve(top);
        if (relax_res.status() == infeasible) break;
        if (relax_res.status() == solved and
            relax_res.optimalValue() <= best) break;
        if (result_is_integral(relax_res)) {
          best = result.obj
        }
      } while();
    }
  }
}
