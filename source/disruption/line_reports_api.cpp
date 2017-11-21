/* Copyright © 2001-2017, Canal TP and/or its affiliates. All rights reserved.
  
This file is part of Navitia,
    the software to build cool stuff with public transport.
 
Hope you'll enjoy and contribute to this project,
    powered by Canal TP (www.canaltp.fr).
Help us simplify mobility and open public transport:
    a non ending quest to the responsive locomotion way of traveling!
  
LICENCE: This program is free software; you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
   
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.
   
You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
  
Stay tuned using
twitter @navitia 
IRC #navitia on freenode
https://groups.google.com/d/forum/navitia
www.navitia.io
*/

#include "line_reports_api.h"

namespace bt = boost::posix_time;
namespace nt = navitia::type;

namespace navitia { namespace disruption {

namespace {

struct LineReport {
    const nt::Line* line;
    std::vector<const nt::Network*> networks;
    std::vector<const nt::Route*> routes;
    std::vector<const nt::StopArea*> stop_areas;
    std::vector<const nt::StopPoint*> stop_points;

    LineReport(const nt::Line* line,
               const std::string& filter,
               const std::vector<std::string>& forbidden_uris,
               const type::Data& d,
               const boost::posix_time::ptime now):
        line(line)
    {
        add_objects(filter, forbidden_uris, d, now, networks);
        add_objects(filter, forbidden_uris, d, now, routes);
        add_objects(filter, forbidden_uris, d, now, stop_areas);
        add_objects(filter, forbidden_uris, d, now, stop_points);
    }

    template<typename T>
    void add_objects(const std::string& filter,
                     const std::vector<std::string>& forbidden_uris,
                     const type::Data& d,
                     const boost::posix_time::ptime now,
                     std::vector<const T*>& objects) {
        std::string new_filter = "line.uri=" + line->uri;
        if (!filter.empty()) { new_filter += " and " + filter; }

        type::Indexes indices;
        try {
            indices = ptref::make_query(nt::get_type_e<T>(), new_filter, forbidden_uris, d);
        } catch (const std::exception&) {}

        for (const auto& idx: indices) {
            const auto* obj = d.pt_data->collection<T>()[idx];
            if (obj->has_publishable_message(now)) {
                objects.push_back(obj);
            }
        }
    }

    bool has_disruption(const boost::posix_time::ptime& current_time) const {
        return line->has_publishable_message(current_time)
            || ! networks.empty()
            || ! routes.empty()
            || ! stop_areas.empty()
            || ! stop_points.empty();
    }

    void to_pb(navitia::PbCreator& pb_creator, const size_t depth) const {
        auto* report = pb_creator.add_line_reports();
        pb_creator.fill(line, report->mutable_line(), depth-1);
        if (line->has_publishable_message(pb_creator.now)) {
            pb_creator.fill(line, report->add_pt_objects(), 0);
        }
        pb_creator.fill(networks, report->mutable_pt_objects(), 0);
        pb_creator.fill(routes, report->mutable_pt_objects(), 0);
        pb_creator.fill(stop_areas, report->mutable_pt_objects(), 0);
        pb_creator.fill(stop_points, report->mutable_pt_objects(), 0);
    }
};

}

void line_reports(navitia::PbCreator& pb_creator,
                  const navitia::type::Data& d,
                  const size_t depth,
                  size_t count,
                  size_t start_page, const std::string& filter,
                  const std::vector<std::string>& forbidden_uris) {
    pb_creator.action_period = bt::time_period(pb_creator.now, bt::seconds(1));
    type::Indexes line_indices;
    try {
        line_indices = ptref::make_query(type::Type_e::Line, filter, forbidden_uris, d);
    } catch(const ptref::parsing_error& parse_error) {
        pb_creator.fill_pb_error(pbnavitia::Error::unable_to_parse, "Unable to parse filter" + parse_error.more);
        return;
    } catch(const ptref::ptref_error& ptref_error) {
        pb_creator.fill_pb_error(pbnavitia::Error::bad_filter, "ptref : "  + ptref_error.more);
        return;
    }
    std::vector<LineReport> line_reports;
    for (auto idx: line_indices) {
        auto line_report = LineReport(d.pt_data->lines[idx], filter, forbidden_uris, d, pb_creator.now);
        if (line_report.has_disruption(pb_creator.now)) {
            line_reports.push_back(line_report);
        }
    }
    const auto total_results = line_reports.size();
    std::vector<LineReport> paged_line_reports = paginate(line_reports, count, start_page);
    std::sort(paged_line_reports.begin(),
    		  paged_line_reports.end(),
              [&](const LineReport& lhs, const LineReport& rhs){
                  return *lhs.line < *rhs.line;});
    for (const auto& line_report: paged_line_reports) {
        line_report.to_pb(pb_creator, depth);
    }
    pb_creator.make_paginate(total_results, start_page, count, pb_creator.line_reports_size());
    if (pb_creator.line_reports_size() == 0) {
        pb_creator.fill_pb_error(pbnavitia::Error::no_solution, pbnavitia::NO_SOLUTION,
                                 "no result for this request");
    }
}

}}// navitia::disruption
