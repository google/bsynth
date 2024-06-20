// Generated by rstantools.  Do not edit by hand.

/*
    bsynth is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    bsynth is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with bsynth.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#ifndef USE_STANC3
#define USE_STANC3
#endif
#include <rstan/rstaninc.hpp>
// Code generated by stanc v2.32.2
#include <stan/model/model_header.hpp>
namespace model_model5_namespace {
using stan::model::model_base_crtp;
using namespace stan::math;
stan::math::profile_map profiles__;
static constexpr std::array<const char*, 66> locations_array__ =
  {" (found before start of program)",
  " (in 'model5', line 51, column 3 to column 26)",
  " (in 'model5', line 52, column 3 to column 22)",
  " (in 'model5', line 62, column 3 to column 22)",
  " (in 'model5', line 63, column 3 to column 28)",
  " (in 'model5', line 66, column 6 to line 67, column 32)",
  " (in 'model5', line 65, column 21 to line 68, column 6)",
  " (in 'model5', line 65, column 6 to line 68, column 6)",
  " (in 'model5', line 70, column 6 to line 71, column 30)",
  " (in 'model5', line 69, column 25 to line 72, column 6)",
  " (in 'model5', line 69, column 5 to line 72, column 6)",
  " (in 'model5', line 64, column 16 to line 73, column 4)",
  " (in 'model5', line 64, column 3 to line 73, column 4)",
  " (in 'model5', line 56, column 6 to column 29)",
  " (in 'model5', line 57, column 6 to line 58, column 38)",
  " (in 'model5', line 55, column 16 to line 59, column 4)",
  " (in 'model5', line 55, column 3 to line 59, column 4)",
  " (in 'model5', line 20, column 3 to column 18)",
  " (in 'model5', line 21, column 3 to column 18)",
  " (in 'model5', line 22, column 15 to column 16)",
  " (in 'model5', line 22, column 10 to column 11)",
  " (in 'model5', line 22, column 3 to column 18)",
  " (in 'model5', line 23, column 3 to column 18)",
  " (in 'model5', line 24, column 10 to column 11)",
  " (in 'model5', line 24, column 12 to column 13)",
  " (in 'model5', line 24, column 3 to column 17)",
  " (in 'model5', line 25, column 3 to column 23)",
  " (in 'model5', line 26, column 10 to column 16)",
  " (in 'model5', line 26, column 17 to column 18)",
  " (in 'model5', line 26, column 3 to column 27)",
  " (in 'model5', line 29, column 10 to column 11)",
  " (in 'model5', line 29, column 13 to column 14)",
  " (in 'model5', line 29, column 3 to column 22)",
  " (in 'model5', line 30, column 10 to column 16)",
  " (in 'model5', line 30, column 18 to column 19)",
  " (in 'model5', line 30, column 3 to column 32)",
  " (in 'model5', line 31, column 19 to column 20)",
  " (in 'model5', line 31, column 10 to column 11)",
  " (in 'model5', line 31, column 3 to column 22)",
  " (in 'model5', line 32, column 15 to column 16)",
  " (in 'model5', line 32, column 3 to column 18)",
  " (in 'model5', line 33, column 12 to column 13)",
  " (in 'model5', line 33, column 2 to column 15)",
  " (in 'model5', line 35, column 12 to column 13)",
  " (in 'model5', line 35, column 5 to column 22)",
  " (in 'model5', line 36, column 11 to column 12)",
  " (in 'model5', line 36, column 4 to column 19)",
  " (in 'model5', line 38, column 8 to column 32)",
  " (in 'model5', line 39, column 8 to column 28)",
  " (in 'model5', line 40, column 8 to column 50)",
  " (in 'model5', line 41, column 8 to column 60)",
  " (in 'model5', line 37, column 19 to line 42, column 8)",
  " (in 'model5', line 37, column 4 to line 42, column 8)",
  " (in 'model5', line 44, column 8 to column 31)",
  " (in 'model5', line 45, column 7 to column 26)",
  " (in 'model5', line 46, column 8 to column 48)",
  " (in 'model5', line 43, column 19 to line 47, column 6)",
  " (in 'model5', line 43, column 5 to line 47, column 6)",
  " (in 'model5', line 34, column 3 to line 48, column 4)",
  " (in 'model5', line 51, column 23 to column 24)",
  " (in 'model5', line 52, column 19 to column 20)",
  " (in 'model5', line 52, column 11 to column 12)",
  " (in 'model5', line 62, column 19 to column 20)",
  " (in 'model5', line 62, column 10 to column 11)",
  " (in 'model5', line 63, column 25 to column 26)",
  " (in 'model5', line 63, column 10 to column 16)"};
#include <stan_meta_header.hpp>
class model_model5 final : public model_base_crtp<model_model5> {
private:
  int N;
  int I;
  std::vector<Eigen::Matrix<double,-1,1>> y;
  int K;
  Eigen::Matrix<double,-1,-1> X_data__;
  int N_pred;
  Eigen::Matrix<double,-1,-1> X_pred_data__;
  Eigen::Matrix<double,-1,-1> X_std_data__;
  Eigen::Matrix<double,-1,-1> X_pred_std_data__;
  std::vector<Eigen::Matrix<double,-1,1>> y_std;
  std::vector<double> mean_y;
  std::vector<double> sd_y;
  Eigen::Map<Eigen::Matrix<double,-1,-1>> X{nullptr, 0, 0};
  Eigen::Map<Eigen::Matrix<double,-1,-1>> X_pred{nullptr, 0, 0};
  Eigen::Map<Eigen::Matrix<double,-1,-1>> X_std{nullptr, 0, 0};
  Eigen::Map<Eigen::Matrix<double,-1,-1>> X_pred_std{nullptr, 0, 0};
public:
  ~model_model5() {}
  model_model5(stan::io::var_context& context__, unsigned int
               random_seed__ = 0, std::ostream* pstream__ = nullptr)
      : model_base_crtp(0) {
    int current_statement__ = 0;
    using local_scalar_t__ = double;
    boost::ecuyer1988 base_rng__ =
      stan::services::util::create_rng(random_seed__, 0);
    // suppress unused var warning
    (void) base_rng__;
    static constexpr const char* function__ =
      "model_model5_namespace::model_model5";
    // suppress unused var warning
    (void) function__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      current_statement__ = 17;
      context__.validate_dims("data initialization", "N", "int",
        std::vector<size_t>{});
      N = std::numeric_limits<int>::min();
      current_statement__ = 17;
      N = context__.vals_i("N")[(1 - 1)];
      current_statement__ = 17;
      stan::math::check_greater_or_equal(function__, "N", N, 1);
      current_statement__ = 18;
      context__.validate_dims("data initialization", "I", "int",
        std::vector<size_t>{});
      I = std::numeric_limits<int>::min();
      current_statement__ = 18;
      I = context__.vals_i("I")[(1 - 1)];
      current_statement__ = 18;
      stan::math::check_greater_or_equal(function__, "I", I, 1);
      current_statement__ = 19;
      stan::math::validate_non_negative_index("y", "I", I);
      current_statement__ = 20;
      stan::math::validate_non_negative_index("y", "N", N);
      current_statement__ = 21;
      context__.validate_dims("data initialization", "y", "double",
        std::vector<size_t>{static_cast<size_t>(I), static_cast<size_t>(N)});
      y = std::vector<Eigen::Matrix<double,-1,1>>(I,
            Eigen::Matrix<double,-1,1>::Constant(N,
              std::numeric_limits<double>::quiet_NaN()));
      {
        std::vector<local_scalar_t__> y_flat__;
        current_statement__ = 21;
        y_flat__ = context__.vals_r("y");
        current_statement__ = 21;
        pos__ = 1;
        current_statement__ = 21;
        for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
          current_statement__ = 21;
          for (int sym2__ = 1; sym2__ <= I; ++sym2__) {
            current_statement__ = 21;
            stan::model::assign(y, y_flat__[(pos__ - 1)],
              "assigning variable y", stan::model::index_uni(sym2__),
              stan::model::index_uni(sym1__));
            current_statement__ = 21;
            pos__ = (pos__ + 1);
          }
        }
      }
      current_statement__ = 22;
      context__.validate_dims("data initialization", "K", "int",
        std::vector<size_t>{});
      K = std::numeric_limits<int>::min();
      current_statement__ = 22;
      K = context__.vals_i("K")[(1 - 1)];
      current_statement__ = 22;
      stan::math::check_greater_or_equal(function__, "K", K, 0);
      current_statement__ = 23;
      stan::math::validate_non_negative_index("X", "N", N);
      current_statement__ = 24;
      stan::math::validate_non_negative_index("X", "K", K);
      current_statement__ = 25;
      context__.validate_dims("data initialization", "X", "double",
        std::vector<size_t>{static_cast<size_t>(N), static_cast<size_t>(K)});
      X_data__ = Eigen::Matrix<double,-1,-1>::Constant(N, K,
                   std::numeric_limits<double>::quiet_NaN());
      new (&X) Eigen::Map<Eigen::Matrix<double,-1,-1>>(X_data__.data(), N, K);
      {
        std::vector<local_scalar_t__> X_flat__;
        current_statement__ = 25;
        X_flat__ = context__.vals_r("X");
        current_statement__ = 25;
        pos__ = 1;
        current_statement__ = 25;
        for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
          current_statement__ = 25;
          for (int sym2__ = 1; sym2__ <= N; ++sym2__) {
            current_statement__ = 25;
            stan::model::assign(X, X_flat__[(pos__ - 1)],
              "assigning variable X", stan::model::index_uni(sym2__),
              stan::model::index_uni(sym1__));
            current_statement__ = 25;
            pos__ = (pos__ + 1);
          }
        }
      }
      current_statement__ = 26;
      context__.validate_dims("data initialization", "N_pred", "int",
        std::vector<size_t>{});
      N_pred = std::numeric_limits<int>::min();
      current_statement__ = 26;
      N_pred = context__.vals_i("N_pred")[(1 - 1)];
      current_statement__ = 26;
      stan::math::check_greater_or_equal(function__, "N_pred", N_pred, 1);
      current_statement__ = 27;
      stan::math::validate_non_negative_index("X_pred", "N_pred", N_pred);
      current_statement__ = 28;
      stan::math::validate_non_negative_index("X_pred", "K", K);
      current_statement__ = 29;
      context__.validate_dims("data initialization", "X_pred", "double",
        std::vector<size_t>{static_cast<size_t>(N_pred),
          static_cast<size_t>(K)});
      X_pred_data__ = Eigen::Matrix<double,-1,-1>::Constant(N_pred, K,
                        std::numeric_limits<double>::quiet_NaN());
      new (&X_pred)
        Eigen::Map<Eigen::Matrix<double,-1,-1>>(X_pred_data__.data(), N_pred,
        K);
      {
        std::vector<local_scalar_t__> X_pred_flat__;
        current_statement__ = 29;
        X_pred_flat__ = context__.vals_r("X_pred");
        current_statement__ = 29;
        pos__ = 1;
        current_statement__ = 29;
        for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
          current_statement__ = 29;
          for (int sym2__ = 1; sym2__ <= N_pred; ++sym2__) {
            current_statement__ = 29;
            stan::model::assign(X_pred, X_pred_flat__[(pos__ - 1)],
              "assigning variable X_pred", stan::model::index_uni(sym2__),
              stan::model::index_uni(sym1__));
            current_statement__ = 29;
            pos__ = (pos__ + 1);
          }
        }
      }
      current_statement__ = 30;
      stan::math::validate_non_negative_index("X_std", "N", N);
      current_statement__ = 31;
      stan::math::validate_non_negative_index("X_std", "K", K);
      current_statement__ = 32;
      X_std_data__ = Eigen::Matrix<double,-1,-1>::Constant(N, K,
                       std::numeric_limits<double>::quiet_NaN());
      new (&X_std)
        Eigen::Map<Eigen::Matrix<double,-1,-1>>(X_std_data__.data(), N, K);
      current_statement__ = 33;
      stan::math::validate_non_negative_index("X_pred_std", "N_pred", N_pred);
      current_statement__ = 34;
      stan::math::validate_non_negative_index("X_pred_std", "K", K);
      current_statement__ = 35;
      X_pred_std_data__ = Eigen::Matrix<double,-1,-1>::Constant(N_pred, K,
                            std::numeric_limits<double>::quiet_NaN());
      new (&X_pred_std)
        Eigen::Map<Eigen::Matrix<double,-1,-1>>(X_pred_std_data__.data(),
        N_pred, K);
      current_statement__ = 36;
      stan::math::validate_non_negative_index("y_std", "I", I);
      current_statement__ = 37;
      stan::math::validate_non_negative_index("y_std", "N", N);
      current_statement__ = 38;
      y_std = std::vector<Eigen::Matrix<double,-1,1>>(I,
                Eigen::Matrix<double,-1,1>::Constant(N,
                  std::numeric_limits<double>::quiet_NaN()));
      current_statement__ = 39;
      stan::math::validate_non_negative_index("mean_y", "I", I);
      current_statement__ = 40;
      mean_y = std::vector<double>(I,
                 std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 41;
      stan::math::validate_non_negative_index("sd_y", "I", I);
      current_statement__ = 42;
      sd_y = std::vector<double>(I, std::numeric_limits<double>::quiet_NaN());
      {
        current_statement__ = 43;
        stan::math::validate_non_negative_index("mean_X", "K", K);
        Eigen::Matrix<double,-1,1> mean_X =
          Eigen::Matrix<double,-1,1>::Constant(K,
            std::numeric_limits<double>::quiet_NaN());
        current_statement__ = 45;
        stan::math::validate_non_negative_index("sd_X", "K", K);
        Eigen::Matrix<double,-1,1> sd_X =
          Eigen::Matrix<double,-1,1>::Constant(K,
            std::numeric_limits<double>::quiet_NaN());
        current_statement__ = 52;
        for (int k = 1; k <= K; ++k) {
          current_statement__ = 47;
          stan::model::assign(mean_X,
            stan::math::mean(
              stan::model::rvalue(X, "X", stan::model::index_omni(),
                stan::model::index_uni(k))), "assigning variable mean_X",
            stan::model::index_uni(k));
          current_statement__ = 48;
          stan::model::assign(sd_X,
            stan::math::sd(
              stan::model::rvalue(X, "X", stan::model::index_omni(),
                stan::model::index_uni(k))), "assigning variable sd_X",
            stan::model::index_uni(k));
          current_statement__ = 49;
          stan::model::assign(X_std,
            stan::math::divide(
              stan::math::subtract(
                stan::model::rvalue(X, "X", stan::model::index_omni(),
                  stan::model::index_uni(k)),
                stan::model::rvalue(mean_X, "mean_X",
                  stan::model::index_uni(k))),
              stan::model::rvalue(sd_X, "sd_X", stan::model::index_uni(k))),
            "assigning variable X_std", stan::model::index_omni(),
            stan::model::index_uni(k));
          current_statement__ = 50;
          stan::model::assign(X_pred_std,
            stan::math::divide(
              stan::math::subtract(
                stan::model::rvalue(X_pred, "X_pred",
                  stan::model::index_omni(), stan::model::index_uni(k)),
                stan::model::rvalue(mean_X, "mean_X",
                  stan::model::index_uni(k))),
              stan::model::rvalue(sd_X, "sd_X", stan::model::index_uni(k))),
            "assigning variable X_pred_std", stan::model::index_omni(),
            stan::model::index_uni(k));
        }
        current_statement__ = 57;
        for (int i = 1; i <= I; ++i) {
          current_statement__ = 53;
          stan::model::assign(mean_y,
            stan::math::mean(
              stan::model::rvalue(y, "y", stan::model::index_uni(i))),
            "assigning variable mean_y", stan::model::index_uni(i));
          current_statement__ = 54;
          stan::model::assign(sd_y,
            stan::math::sd(
              stan::model::rvalue(y, "y", stan::model::index_uni(i))),
            "assigning variable sd_y", stan::model::index_uni(i));
          current_statement__ = 55;
          stan::model::assign(y_std,
            stan::math::divide(
              stan::math::subtract(
                stan::model::rvalue(y, "y", stan::model::index_uni(i)),
                stan::model::rvalue(mean_y, "mean_y",
                  stan::model::index_uni(i))),
              stan::model::rvalue(sd_y, "sd_y", stan::model::index_uni(i))),
            "assigning variable y_std", stan::model::index_uni(i));
        }
      }
      current_statement__ = 59;
      stan::math::validate_non_negative_index("sigma", "I", I);
      current_statement__ = 60;
      stan::math::validate_non_negative_index("beta", "I", I);
      current_statement__ = 61;
      stan::math::validate_positive_index("beta", "K", K);
      current_statement__ = 62;
      stan::math::validate_non_negative_index("y_sim", "I", I);
      current_statement__ = 63;
      stan::math::validate_non_negative_index("y_sim", "N", N);
      current_statement__ = 64;
      stan::math::validate_non_negative_index("y_pred", "I", I);
      current_statement__ = 65;
      stan::math::validate_non_negative_index("y_pred", "N_pred", N_pred);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    num_params_r__ = I + (I * (K - 1));
  }
  inline std::string model_name() const final {
    return "model_model5";
  }
  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 v2.32.2",
             "stancflags = --allow-undefined"};
  }
  template <bool propto__, bool jacobian__, typename VecR, typename VecI,
            stan::require_vector_like_t<VecR>* = nullptr,
            stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline stan::scalar_type_t<VecR>
  log_prob_impl(VecR& params_r__, VecI& params_i__, std::ostream*
                pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    static constexpr const char* function__ =
      "model_model5_namespace::log_prob";
    // suppress unused var warning
    (void) function__;
    try {
      std::vector<local_scalar_t__> sigma =
        std::vector<local_scalar_t__>(I, DUMMY_VAR__);
      current_statement__ = 1;
      sigma = in__.template read_constrain_lb<std::vector<local_scalar_t__>,
                jacobian__>(0, lp__, I);
      std::vector<Eigen::Matrix<local_scalar_t__,-1,1>> beta =
        std::vector<Eigen::Matrix<local_scalar_t__,-1,1>>(I,
          Eigen::Matrix<local_scalar_t__,-1,1>::Constant(K, DUMMY_VAR__));
      current_statement__ = 2;
      beta = in__.template read_constrain_simplex<
               std::vector<Eigen::Matrix<local_scalar_t__,-1,1>>,
               jacobian__>(lp__, I, K);
      {
        current_statement__ = 16;
        for (int i = 1; i <= I; ++i) {
          current_statement__ = 13;
          lp_accum__.add(stan::math::normal_lpdf<propto__>(
                           stan::model::rvalue(sigma, "sigma",
                             stan::model::index_uni(i)), 0, 1));
          current_statement__ = 14;
          lp_accum__.add(stan::math::normal_lpdf<false>(
                           stan::model::rvalue(y_std, "y_std",
                             stan::model::index_uni(i)),
                           stan::math::multiply(X_std,
                             stan::model::rvalue(beta, "beta",
                               stan::model::index_uni(i))),
                           stan::model::rvalue(sigma, "sigma",
                             stan::model::index_uni(i))));
        }
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
  }
  template <typename RNG, typename VecR, typename VecI, typename VecVar,
            stan::require_vector_like_vt<std::is_floating_point,
            VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral,
            VecI>* = nullptr, stan::require_vector_vt<std::is_floating_point,
            VecVar>* = nullptr>
  inline void
  write_array_impl(RNG& base_rng__, VecR& params_r__, VecI& params_i__,
                   VecVar& vars__, const bool
                   emit_transformed_parameters__ = true, const bool
                   emit_generated_quantities__ = true, std::ostream*
                   pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    static constexpr bool propto__ = true;
    // suppress unused var warning
    (void) propto__;
    double lp__ = 0.0;
    // suppress unused var warning
    (void) lp__;
    int current_statement__ = 0;
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    constexpr bool jacobian__ = false;
    static constexpr const char* function__ =
      "model_model5_namespace::write_array";
    // suppress unused var warning
    (void) function__;
    try {
      std::vector<double> sigma =
        std::vector<double>(I, std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 1;
      sigma = in__.template read_constrain_lb<std::vector<local_scalar_t__>,
                jacobian__>(0, lp__, I);
      std::vector<Eigen::Matrix<double,-1,1>> beta =
        std::vector<Eigen::Matrix<double,-1,1>>(I,
          Eigen::Matrix<double,-1,1>::Constant(K,
            std::numeric_limits<double>::quiet_NaN()));
      current_statement__ = 2;
      beta = in__.template read_constrain_simplex<
               std::vector<Eigen::Matrix<local_scalar_t__,-1,1>>,
               jacobian__>(lp__, I, K);
      out__.write(sigma);
      for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
        for (int sym2__ = 1; sym2__ <= I; ++sym2__) {
          out__.write(beta[(sym2__ - 1)][(sym1__ - 1)]);
        }
      }
      if (stan::math::logical_negation(
            (stan::math::primitive_value(emit_transformed_parameters__) ||
            stan::math::primitive_value(emit_generated_quantities__)))) {
        return ;
      }
      if (stan::math::logical_negation(emit_generated_quantities__)) {
        return ;
      }
      std::vector<Eigen::Matrix<double,-1,1>> y_sim =
        std::vector<Eigen::Matrix<double,-1,1>>(I,
          Eigen::Matrix<double,-1,1>::Constant(N,
            std::numeric_limits<double>::quiet_NaN()));
      std::vector<Eigen::Matrix<double,-1,1>> y_pred =
        std::vector<Eigen::Matrix<double,-1,1>>(I,
          Eigen::Matrix<double,-1,1>::Constant(N_pred,
            std::numeric_limits<double>::quiet_NaN()));
      current_statement__ = 12;
      for (int i = 1; i <= I; ++i) {
        current_statement__ = 7;
        for (int n = 1; n <= N; ++n) {
          current_statement__ = 5;
          stan::model::assign(y_sim,
            ((stan::math::normal_rng(
                stan::math::multiply(
                  stan::model::rvalue(X_std, "X_std",
                    stan::model::index_uni(n), stan::model::index_omni()),
                  stan::model::rvalue(beta, "beta", stan::model::index_uni(i))),
                stan::model::rvalue(sigma, "sigma", stan::model::index_uni(i)),
                base_rng__) *
            stan::model::rvalue(sd_y, "sd_y", stan::model::index_uni(i))) +
            stan::model::rvalue(mean_y, "mean_y", stan::model::index_uni(i))),
            "assigning variable y_sim", stan::model::index_uni(i),
            stan::model::index_uni(n));
        }
        current_statement__ = 10;
        for (int j = 1; j <= N_pred; ++j) {
          current_statement__ = 8;
          stan::model::assign(y_pred,
            ((stan::math::normal_rng(
                stan::math::multiply(
                  stan::model::rvalue(X_pred_std, "X_pred_std",
                    stan::model::index_uni(j), stan::model::index_omni()),
                  stan::model::rvalue(beta, "beta", stan::model::index_uni(i))),
                stan::model::rvalue(sigma, "sigma", stan::model::index_uni(i)),
                base_rng__) *
            stan::model::rvalue(sd_y, "sd_y", stan::model::index_uni(i))) +
            stan::model::rvalue(mean_y, "mean_y", stan::model::index_uni(i))),
            "assigning variable y_pred", stan::model::index_uni(i),
            stan::model::index_uni(j));
        }
      }
      for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
        for (int sym2__ = 1; sym2__ <= I; ++sym2__) {
          out__.write(y_sim[(sym2__ - 1)][(sym1__ - 1)]);
        }
      }
      for (int sym1__ = 1; sym1__ <= N_pred; ++sym1__) {
        for (int sym2__ = 1; sym2__ <= I; ++sym2__) {
          out__.write(y_pred[(sym2__ - 1)][(sym1__ - 1)]);
        }
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  template <typename VecVar, typename VecI,
            stan::require_vector_t<VecVar>* = nullptr,
            stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline void
  unconstrain_array_impl(const VecVar& params_r__, const VecI& params_i__,
                         VecVar& vars__, std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      std::vector<local_scalar_t__> sigma =
        std::vector<local_scalar_t__>(I, DUMMY_VAR__);
      current_statement__ = 1;
      stan::model::assign(sigma, in__.read<std::vector<local_scalar_t__>>(I),
        "assigning variable sigma");
      out__.write_free_lb(0, sigma);
      std::vector<Eigen::Matrix<local_scalar_t__,-1,1>> beta =
        std::vector<Eigen::Matrix<local_scalar_t__,-1,1>>(I,
          Eigen::Matrix<local_scalar_t__,-1,1>::Constant(K, DUMMY_VAR__));
      current_statement__ = 2;
      for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
        current_statement__ = 2;
        for (int sym2__ = 1; sym2__ <= I; ++sym2__) {
          current_statement__ = 2;
          stan::model::assign(beta, in__.read<local_scalar_t__>(),
            "assigning variable beta", stan::model::index_uni(sym2__),
            stan::model::index_uni(sym1__));
        }
      }
      out__.write_free_simplex(beta);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  template <typename VecVar, stan::require_vector_t<VecVar>* = nullptr>
  inline void
  transform_inits_impl(const stan::io::var_context& context__, VecVar&
                       vars__, std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      current_statement__ = 1;
      context__.validate_dims("parameter initialization", "sigma", "double",
        std::vector<size_t>{static_cast<size_t>(I)});
      current_statement__ = 2;
      context__.validate_dims("parameter initialization", "beta", "double",
        std::vector<size_t>{static_cast<size_t>(I), static_cast<size_t>(K)});
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      std::vector<local_scalar_t__> sigma =
        std::vector<local_scalar_t__>(I, DUMMY_VAR__);
      current_statement__ = 1;
      sigma = context__.vals_r("sigma");
      out__.write_free_lb(0, sigma);
      std::vector<Eigen::Matrix<local_scalar_t__,-1,1>> beta =
        std::vector<Eigen::Matrix<local_scalar_t__,-1,1>>(I,
          Eigen::Matrix<local_scalar_t__,-1,1>::Constant(K, DUMMY_VAR__));
      {
        std::vector<local_scalar_t__> beta_flat__;
        current_statement__ = 2;
        beta_flat__ = context__.vals_r("beta");
        current_statement__ = 2;
        pos__ = 1;
        current_statement__ = 2;
        for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
          current_statement__ = 2;
          for (int sym2__ = 1; sym2__ <= I; ++sym2__) {
            current_statement__ = 2;
            stan::model::assign(beta, beta_flat__[(pos__ - 1)],
              "assigning variable beta", stan::model::index_uni(sym2__),
              stan::model::index_uni(sym1__));
            current_statement__ = 2;
            pos__ = (pos__ + 1);
          }
        }
      }
      out__.write_free_simplex(beta);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  inline void
  get_param_names(std::vector<std::string>& names__, const bool
                  emit_transformed_parameters__ = true, const bool
                  emit_generated_quantities__ = true) const {
    names__ = std::vector<std::string>{"sigma", "beta"};
    if (emit_transformed_parameters__) {}
    if (emit_generated_quantities__) {
      std::vector<std::string> temp{"y_sim", "y_pred"};
      names__.reserve(names__.size() + temp.size());
      names__.insert(names__.end(), temp.begin(), temp.end());
    }
  }
  inline void
  get_dims(std::vector<std::vector<size_t>>& dimss__, const bool
           emit_transformed_parameters__ = true, const bool
           emit_generated_quantities__ = true) const {
    dimss__ = std::vector<std::vector<size_t>>{std::vector<size_t>{static_cast<
                                                                    size_t>(I)},
                std::vector<size_t>{static_cast<size_t>(I),
                  static_cast<size_t>(K)}};
    if (emit_transformed_parameters__) {}
    if (emit_generated_quantities__) {
      std::vector<std::vector<size_t>>
        temp{std::vector<size_t>{static_cast<size_t>(I),
               static_cast<size_t>(N)},
             std::vector<size_t>{static_cast<size_t>(I),
               static_cast<size_t>(N_pred)}};
      dimss__.reserve(dimss__.size() + temp.size());
      dimss__.insert(dimss__.end(), temp.begin(), temp.end());
    }
  }
  inline void
  constrained_param_names(std::vector<std::string>& param_names__, bool
                          emit_transformed_parameters__ = true, bool
                          emit_generated_quantities__ = true) const final {
    for (int sym1__ = 1; sym1__ <= I; ++sym1__) {
      param_names__.emplace_back(std::string() + "sigma" + '.' +
        std::to_string(sym1__));
    }
    for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
      for (int sym2__ = 1; sym2__ <= I; ++sym2__) {
        param_names__.emplace_back(std::string() + "beta" + '.' +
          std::to_string(sym2__) + '.' + std::to_string(sym1__));
      }
    }
    if (emit_transformed_parameters__) {}
    if (emit_generated_quantities__) {
      for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
        for (int sym2__ = 1; sym2__ <= I; ++sym2__) {
          param_names__.emplace_back(std::string() + "y_sim" + '.' +
            std::to_string(sym2__) + '.' + std::to_string(sym1__));
        }
      }
      for (int sym1__ = 1; sym1__ <= N_pred; ++sym1__) {
        for (int sym2__ = 1; sym2__ <= I; ++sym2__) {
          param_names__.emplace_back(std::string() + "y_pred" + '.' +
            std::to_string(sym2__) + '.' + std::to_string(sym1__));
        }
      }
    }
  }
  inline void
  unconstrained_param_names(std::vector<std::string>& param_names__, bool
                            emit_transformed_parameters__ = true, bool
                            emit_generated_quantities__ = true) const final {
    for (int sym1__ = 1; sym1__ <= I; ++sym1__) {
      param_names__.emplace_back(std::string() + "sigma" + '.' +
        std::to_string(sym1__));
    }
    for (int sym1__ = 1; sym1__ <= (K - 1); ++sym1__) {
      for (int sym2__ = 1; sym2__ <= I; ++sym2__) {
        param_names__.emplace_back(std::string() + "beta" + '.' +
          std::to_string(sym2__) + '.' + std::to_string(sym1__));
      }
    }
    if (emit_transformed_parameters__) {}
    if (emit_generated_quantities__) {
      for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
        for (int sym2__ = 1; sym2__ <= I; ++sym2__) {
          param_names__.emplace_back(std::string() + "y_sim" + '.' +
            std::to_string(sym2__) + '.' + std::to_string(sym1__));
        }
      }
      for (int sym1__ = 1; sym1__ <= N_pred; ++sym1__) {
        for (int sym2__ = 1; sym2__ <= I; ++sym2__) {
          param_names__.emplace_back(std::string() + "y_pred" + '.' +
            std::to_string(sym2__) + '.' + std::to_string(sym1__));
        }
      }
    }
  }
  inline std::string get_constrained_sizedtypes() const {
    return std::string("[{\"name\":\"sigma\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(I) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"parameters\"},{\"name\":\"beta\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(I) + ",\"element_type\":{\"name\":\"vector\",\"length\":" + std::to_string(K) + "}},\"block\":\"parameters\"},{\"name\":\"y_sim\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(I) + ",\"element_type\":{\"name\":\"vector\",\"length\":" + std::to_string(N) + "}},\"block\":\"generated_quantities\"},{\"name\":\"y_pred\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(I) + ",\"element_type\":{\"name\":\"vector\",\"length\":" + std::to_string(N_pred) + "}},\"block\":\"generated_quantities\"}]");
  }
  inline std::string get_unconstrained_sizedtypes() const {
    return std::string("[{\"name\":\"sigma\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(I) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"parameters\"},{\"name\":\"beta\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(I) + ",\"element_type\":{\"name\":\"vector\",\"length\":" + std::to_string((K -1)) + "}},\"block\":\"parameters\"},{\"name\":\"y_sim\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(I) + ",\"element_type\":{\"name\":\"vector\",\"length\":" + std::to_string(N) + "}},\"block\":\"generated_quantities\"},{\"name\":\"y_pred\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(I) + ",\"element_type\":{\"name\":\"vector\",\"length\":" + std::to_string(N_pred) + "}},\"block\":\"generated_quantities\"}]");
  }
  // Begin method overload boilerplate
  template <typename RNG> inline void
  write_array(RNG& base_rng, Eigen::Matrix<double,-1,1>& params_r,
              Eigen::Matrix<double,-1,1>& vars, const bool
              emit_transformed_parameters = true, const bool
              emit_generated_quantities = true, std::ostream*
              pstream = nullptr) const {
    const size_t num_params__ = (I + (I * K));
    const size_t num_transformed = emit_transformed_parameters * (0);
    const size_t num_gen_quantities = emit_generated_quantities * (((I * N) +
      (I * N_pred)));
    const size_t num_to_write = num_params__ + num_transformed +
      num_gen_quantities;
    std::vector<int> params_i;
    vars = Eigen::Matrix<double,-1,1>::Constant(num_to_write,
             std::numeric_limits<double>::quiet_NaN());
    write_array_impl(base_rng, params_r, params_i, vars,
      emit_transformed_parameters, emit_generated_quantities, pstream);
  }
  template <typename RNG> inline void
  write_array(RNG& base_rng, std::vector<double>& params_r, std::vector<int>&
              params_i, std::vector<double>& vars, bool
              emit_transformed_parameters = true, bool
              emit_generated_quantities = true, std::ostream*
              pstream = nullptr) const {
    const size_t num_params__ = (I + (I * K));
    const size_t num_transformed = emit_transformed_parameters * (0);
    const size_t num_gen_quantities = emit_generated_quantities * (((I * N) +
      (I * N_pred)));
    const size_t num_to_write = num_params__ + num_transformed +
      num_gen_quantities;
    vars = std::vector<double>(num_to_write,
             std::numeric_limits<double>::quiet_NaN());
    write_array_impl(base_rng, params_r, params_i, vars,
      emit_transformed_parameters, emit_generated_quantities, pstream);
  }
  template <bool propto__, bool jacobian__, typename T_> inline T_
  log_prob(Eigen::Matrix<T_,-1,1>& params_r, std::ostream* pstream = nullptr) const {
    Eigen::Matrix<int,-1,1> params_i;
    return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
  }
  template <bool propto__, bool jacobian__, typename T_> inline T_
  log_prob(std::vector<T_>& params_r, std::vector<int>& params_i,
           std::ostream* pstream = nullptr) const {
    return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
  }
  inline void
  transform_inits(const stan::io::var_context& context,
                  Eigen::Matrix<double,-1,1>& params_r, std::ostream*
                  pstream = nullptr) const final {
    std::vector<double> params_r_vec(params_r.size());
    std::vector<int> params_i;
    transform_inits(context, params_i, params_r_vec, pstream);
    params_r = Eigen::Map<Eigen::Matrix<double,-1,1>>(params_r_vec.data(),
                 params_r_vec.size());
  }
  inline void
  transform_inits(const stan::io::var_context& context, std::vector<int>&
                  params_i, std::vector<double>& vars, std::ostream*
                  pstream__ = nullptr) const {
    vars.resize(num_params_r__);
    transform_inits_impl(context, vars, pstream__);
  }
  inline void
  unconstrain_array(const std::vector<double>& params_constrained,
                    std::vector<double>& params_unconstrained, std::ostream*
                    pstream = nullptr) const {
    const std::vector<int> params_i;
    params_unconstrained = std::vector<double>(num_params_r__,
                             std::numeric_limits<double>::quiet_NaN());
    unconstrain_array_impl(params_constrained, params_i,
      params_unconstrained, pstream);
  }
  inline void
  unconstrain_array(const Eigen::Matrix<double,-1,1>& params_constrained,
                    Eigen::Matrix<double,-1,1>& params_unconstrained,
                    std::ostream* pstream = nullptr) const {
    const std::vector<int> params_i;
    params_unconstrained = Eigen::Matrix<double,-1,1>::Constant(num_params_r__,
                             std::numeric_limits<double>::quiet_NaN());
    unconstrain_array_impl(params_constrained, params_i,
      params_unconstrained, pstream);
  }
};
}
using stan_model = model_model5_namespace::model_model5;
#ifndef USING_R
// Boilerplate
stan::model::model_base&
new_model(stan::io::var_context& data_context, unsigned int seed,
          std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
stan::math::profile_map& get_stan_profile_data() {
  return model_model5_namespace::profiles__;
}
#endif
#endif
