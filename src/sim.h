//
//  Copyright (c) 2022 - 2023, Thijs Janzen
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <memory>
#include <vector>
#include <array>
#include <random>
#include <tuple>
#include <string>
#include <map>
#include <algorithm>


using num_mat     = std::vector< std::vector<double >>;
using num_mat_mat = std::vector<num_mat>;
using vec_dist    = std::vector< std::discrete_distribution<> >;

enum event_type {speciation, extinction, max_num};

enum finish_type {done, extinct, overshoot, conditioning, not_run_yet,
                  max_types};

struct ltab_species {
  ltab_species(double brts, int parent, int ID) :
   bt_(brts), parent_id_(parent), self_id_(ID), extinct_time_(-1) {
  }

  double get_id() const {
    return(self_id_);
  }

  double get_parent() const  {
    return(parent_id_);
  }

  void set_death(double d) {
    extinct_time_ = d;
  }

  bool is_dead() const {
    if (extinct_time_ > -1) return true;
    return false;
  }

  std::array<double, 5> get_data() const {
    return std::array<double, 5>{bt_, parent_id_, self_id_, extinct_time_};
  }

private:
  const double bt_;
  const double parent_id_;
  const double self_id_;
  double extinct_time_;
};


struct secsse_sim {
  std::mt19937_64 rndgen_;

  std::vector< ltab_species > L;

  std::array<int, 2> track_crowns;
  const std::array<double, 2> parameters;
  std::array<double, 2> rates;


  // external data:
  const double max_t;
  const size_t max_spec;
  const bool non_extinction;
  const bool max_spec_extant;
  const bool crown_start;

  finish_type run_info;
  double t;

  secsse_sim(const std::array<double, 2>& params,
             double total_time,
             size_t max_s,
             bool max_s_e,
             const bool& ne,
             int seed,
             bool start_at_crown) :
    parameters(params),
    max_t(total_time),
    max_spec(max_s),
    non_extinction(ne),
    max_spec_extant(max_s_e),
    crown_start(start_at_crown),
    run_info(not_run_yet),
    t(0.0) {
    // randomize randomizer
    rndgen_.seed((seed < 0) ? std::random_device {}() : seed);
    track_crowns = {0, 0};
    rates = {0.0, 0.0};
  }

  void run() {
    t = 0.0;

    // randomly draw initial trait
    run_info = not_run_yet;

    L.clear();

    if (crown_start) {
      L.push_back(ltab_species(0.0,  0, -1));
      L.push_back(ltab_species(0.0, -1,  2));
    } else {
      L.push_back(ltab_species(0.0,  0, -1));
      track_crowns = {1, 0};
      evolve_until_crown();
      if (t > max_t) {
        run_info = done;
        return;
      }
      if (track_crowns[0] + track_crowns[1] < 1) {
        run_info = extinct;
        return;
      }
    }

    track_crowns = {1, 1};

    while (true) {
      update_rates();
      auto dt = draw_dt();
      t += dt;

      if (t > max_t)  {
        run_info = done; break;
      }

      event_type event = draw_event();
      apply_event(event);

      if (track_crowns[0] < 1 || track_crowns[1] < 1) {
        run_info = extinct;
        break;
      }

      if (max_spec_extant) {
        if (track_crowns[0] + track_crowns[1] >= max_spec) {
          run_info = overshoot; break;
        }
      } else {
        if (L.size() >= max_spec) {
          run_info = overshoot; break;
        }
      }
    }
  }

  void apply_event(const event_type event) {
    switch (event) {
    case speciation: {
      event_speciation();
      break;
    }
    case extinction: {
      event_extinction();
      break;
    }
    default: break;
    }
    return;
  }

  void event_extinction() {
    size_t dying = sample_from_pop();

    if (L[dying].get_id() < 0) {
      track_crowns[0]--;
    } else {
      track_crowns[1]--;
    }

    L[dying].set_death(t);
  }

  void event_speciation() {
    size_t mother = sample_from_pop();

    int mother_id = static_cast<int>(L[mother].get_id());
    if (mother_id == 0) {
      throw "impossible mother";
    }

    int new_id = static_cast<int>(L.size()) + 1;
    if (mother_id < 0) {
      track_crowns[0]++;
      new_id *= -1;
    } else {
      track_crowns[1]++;
    }

    L.emplace_back(ltab_species(t, mother_id, new_id));
  }


  void update_rates() {
    size_t num_species = track_crowns[0] + track_crowns[1];
    rates[event_type::speciation] = num_species                * parameters[event_type::speciation];
    rates[event_type::extinction] = num_species                * parameters[event_type::extinction];
  }

  event_type draw_event() {
    double total_rate = rates[event_type::extinction] +
                        rates[event_type::speciation];

    std::uniform_real_distribution<double> unif_dist(0.0, total_rate);
    double r = unif_dist(rndgen_);

    // ordering of rates is:
    // {shift, speciation, extinction, max_num};
    if (r < rates[event_type::speciation]) return speciation;

    return extinction;
  }

  double draw_dt() {
    double total_rate = rates[event_type::extinction] +
                        rates[event_type::speciation];
    std::exponential_distribution<double> exp_dist(total_rate);
    return exp_dist(rndgen_);
  }

  size_t num_species() {
    if (max_spec_extant) {
      return track_crowns[0] + track_crowns[1];
    } else {
      return L.size();
    }
  }

  size_t ltable_size() {
    return L.size();
  }

  num_mat extract_ltable() {
    num_mat extracted_ltable(L.size(), std::vector<double>(4));
    for (size_t i = 0; i < L.size(); ++i) {
      auto temp = L[i].get_data();
      std::vector<double> row(temp.begin(), temp.end());
      extracted_ltable[i] = row;
    }
    return extracted_ltable;
  }

  void update_tree_size_hist(int* val) {
    if (run_info == extinct) {
      *val = 0;
      return;
    }

    if (max_spec_extant) {
      *val = (track_crowns[0] + track_crowns[1]);
    } else {
      *val = L.size();
    }
    return;
  }

  size_t sample_from_pop() {
    std::uniform_int_distribution<size_t> d(0, L.size() - 1);
    auto index = d(rndgen_);
    while(L[index].is_dead()) {
      index = d(rndgen_);
    }
    return index;
  }

  void evolve_until_crown() {
    while (L.size() < 2) {
      update_rates();
      double dt = draw_dt();
      t += dt;
      if (t > max_t) {
        run_info = done;
        return;
      }
      event_type event = draw_event();
      switch (event) {
      case extinction: {
        event_extinction();
        break;
      }
      case speciation: {
        size_t mother = 0;

        L.emplace_back(t, L[mother].get_id(), 2);
        break;
      }
      default:
        break;
      }

      if (track_crowns[0] + track_crowns[1] < 1) break;
    }
  }
};
