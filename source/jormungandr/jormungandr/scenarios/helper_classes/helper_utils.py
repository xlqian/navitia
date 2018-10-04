# Copyright (c) 2001-2017, Canal TP and/or its affiliates. All rights reserved.
#
# This file is part of Navitia,
#     the software to build cool stuff with public transport.
#
# Hope you'll enjoy and contribute to this project,
#     powered by Canal TP (www.canaltp.fr).
# Help us simplify mobility and open public transport:
#     a non ending quest to the responsive locomotion way of traveling!
#
# LICENCE: This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# Stay tuned using
# twitter @navitia
# IRC #navitia on freenode
# https://groups.google.com/d/forum/navitia
# www.navitia.io

from __future__ import absolute_import
import math
from jormungandr.street_network.street_network import StreetNetworkPathType
from jormungandr.utils import PeriodExtremity, SectionSorter, get_pt_object_coord, generate_id
from navitiacommon import response_pb2
from .helper_exceptions import *
import copy
import logging
import six

MODE_TO_PB_MODE = {
    'walking': response_pb2.Walking,
    'bike': response_pb2.Bike,
    'bss': response_pb2.Bss,
    'car': response_pb2.Car,
    'ridesharing': response_pb2.Ridesharing,
}


def _create_crowfly(pt_journey, crowfly_origin, crowfly_destination, begin, end, mode):
    section = response_pb2.Section()
    section.type = response_pb2.CROW_FLY
    section.origin.CopyFrom(crowfly_origin)
    section.destination.CopyFrom(crowfly_destination)
    section.duration = end - begin
    pt_journey.durations.total += section.duration
    pt_journey.duration += section.duration
    section.begin_date_time = begin
    section.end_date_time = end
    if section.duration > 0:
        section.street_network.mode = MODE_TO_PB_MODE.get(mode)
    # mode is always walking for a teleportation crow_fly
    else:
        section.street_network.mode = response_pb2.Walking

    # Calculate section length
    from_coord = get_pt_object_coord(section.origin)
    to_coord = get_pt_object_coord(section.destination)
    section.length = int(crowfly_distance_between(from_coord, to_coord))

    # The section "distances" and "durations" in the response needs to be updated according to the mode.
    # only if it isn't a 'free' crow_fly
    if section.duration > 0:
        setattr(pt_journey.distances, mode, (getattr(pt_journey.distances, mode) + section.length))
        setattr(pt_journey.durations, mode, (getattr(pt_journey.durations, mode) + section.duration))

    section.id = six.text_type(generate_id())
    return section


def _is_crowfly_needed(uri, fallback_durations, crowfly_sps, fallback_direct_path):
    # Unable to project
    f = fallback_durations.get(uri, None)
    is_unknown_projection = f.status == response_pb2.unknown if f else False

    is_crowfly_sp = uri in crowfly_sps

    # At this point, theoretically, fallback_dp should be found since the isochrone has already given a
    # valid value BUT, in some cases(due to the bad projection, etc), fallback_dp may not exist even
    # though the isochrone said "OK". In these cases, we create a crowfly section
    is_empty_fallback_dp = not fallback_direct_path or not fallback_direct_path.journeys
    return is_unknown_projection or is_crowfly_sp or is_empty_fallback_dp


def _align_fallback_direct_path_datetime(fallback_direct_path, fallback_extremity):
    """
    :param fallback_extremity: is a PeriodExtremity (a datetime and its meaning on the fallback period)
    in this function, we retrieve from streetnetwork_path_pool the direct path regarding to the given
    parameters(mode, orig_uri, etc...) then we recompute the datetimes of the found direct path,
    since the request datetime might no longer be the same (but we consider the same fallback duration).
    """

    # align datetimes to requested ones (as we consider fallback duration are the same no matter when)
    fallback_copy = copy.deepcopy(fallback_direct_path)
    journey = fallback_copy.journeys[0]
    datetime, represents_start_fallback = fallback_extremity
    if represents_start_fallback:
        journey.departure_date_time = datetime
        journey.arrival_date_time = datetime + journey.duration
    else:
        journey.departure_date_time = datetime - journey.duration
        journey.arrival_date_time = datetime

    delta = journey.departure_date_time - journey.sections[0].begin_date_time
    if delta != 0:
        for s in journey.sections:
            s.begin_date_time += delta
            s.end_date_time += delta
    return fallback_copy


def _rename_fallback_sections_ids(sections):
    for s in sections:
        s.id = six.text_type(generate_id())


def _extend_pt_sections_with_fallback_sections(pt_journey, dp_journey):
    if getattr(dp_journey, 'journeys', []) and hasattr(dp_journey.journeys[0], 'sections'):
        _rename_fallback_sections_ids(dp_journey.journeys[0].sections)
        pt_journey.sections.extend(dp_journey.journeys[0].sections)


def _update_journey(pt_journey, aligned_fallback):
    pt_journey.durations.walking += aligned_fallback.journeys[0].durations.walking
    pt_journey.durations.bike += aligned_fallback.journeys[0].durations.bike
    pt_journey.durations.car += aligned_fallback.journeys[0].durations.car

    pt_journey.distances.walking += aligned_fallback.journeys[0].distances.walking
    pt_journey.distances.bike += aligned_fallback.journeys[0].distances.bike
    pt_journey.distances.car += aligned_fallback.journeys[0].distances.car


def _replace_crowfly_with_streetnetwork(pt_journey, fallback_dp, fallback_period_extremity):
    aligned_fallback = _align_fallback_direct_path_datetime(fallback_dp, fallback_period_extremity)

    _update_journey(pt_journey, aligned_fallback)

    if fallback_period_extremity.represents_start:
        section_to_update = pt_journey.sections[-1]
    else:
        section_to_update = pt_journey.sections[0]

    assert section_to_update.type == response_pb2.CROW_FLY
    section_to_update.CopyFrom(aligned_fallback.journeys[0].sections[0])
    assert section_to_update.type == response_pb2.STREET_NETWORK

    ## update the 'id' wich isn't set
    _rename_fallback_sections_ids([section_to_update])


def _extend_journey(pt_journey, fallback_dp, fallback_period_extremity):
    """
    :param fallback_period_extremity: is a PeriodExtremity (a datetime and it's meaning on the fallback period)
    """
    aligned_fallback = _align_fallback_direct_path_datetime(fallback_dp, fallback_period_extremity)

    pt_journey.duration += aligned_fallback.journeys[0].duration
    pt_journey.durations.total = pt_journey.duration
    pt_journey.durations.walking += aligned_fallback.journeys[0].durations.walking
    pt_journey.durations.bike += aligned_fallback.journeys[0].durations.bike
    pt_journey.durations.car += aligned_fallback.journeys[0].durations.car

    pt_journey.distances.walking += aligned_fallback.journeys[0].distances.walking
    pt_journey.distances.bike += aligned_fallback.journeys[0].distances.bike
    pt_journey.distances.car += aligned_fallback.journeys[0].distances.car

    # For start fallback section copy pt_section.origin to last fallback_section.destination
    # where as for end fallback section copy last pt_section.destination to fallback_section.origin
    if fallback_period_extremity.represents_start:
        aligned_fallback.journeys[0].sections[0].origin.CopyFrom(pt_journey.sections[-1].destination)
    else:
        aligned_fallback.journeys[0].sections[-1].destination.CopyFrom(pt_journey.sections[0].origin)

    _extend_pt_sections_with_fallback_sections(pt_journey, aligned_fallback)


def _get_fallback_logic(fallback_type):
    if fallback_type == StreetNetworkPathType.BEGINNING_FALLBACK:
        return BeginningFallback()

    if fallback_type == StreetNetworkPathType.ENDING_FALLBACK:
        return EndingFallback()

    raise EntryPointException(
        "Unknown fallback type : {}".format(fallback_type), response_pb2.Error.unknown_object
    )


class BeginningFallback:
    def get_section(self, pt_journey):
        for s in pt_journey.sections:
            if s.type == response_pb2.CROW_FLY:
                continue
            return s
        return None

    def get_pt_boundaries(self, pt_journey):
        return self.get_section(pt_journey).origin

    def get_pt_bound_datetime(self, pt_journey):
        return pt_journey.departure_date_time

    def represent_start(self):
        return False

    def get_section_datetime(self, pt_journey):
        return self.get_section(pt_journey).begin_date_time

    def get_crowfly_datetime(self, pt_datetime, fallback_duration):
        return pt_datetime - fallback_duration

    def route_params(self, origin, dest):
        return origin, dest

    def set_fallback_datetime(self, pt_journey):
        pt_journey.departure_date_time = self.get_section_datetime(pt_journey)


class EndingFallback:
    def get_section(self, pt_journey):
        for s in reversed(pt_journey.sections):
            if s.type == response_pb2.CROW_FLY:
                continue
            return s
        return None

    def get_pt_boundaries(self, pt_journey):
        return self.get_section(pt_journey).destination

    def get_pt_bound_datetime(self, pt_journey):
        return pt_journey.arrival_date_time

    def represent_start(self):
        return True

    def get_section_datetime(self, pt_journey):
        return self.get_section(pt_journey).end_date_time

    def get_crowfly_datetime(self, pt_datetime, fallback_duration):
        return pt_datetime + fallback_duration

    def route_params(self, origin, dest):
        return dest, origin

    def set_fallback_datetime(self, pt_journey):
        pt_journey.arrival_date_time = self.get_section_datetime(pt_journey)


def _build_crowflies(pt_journeys, orig, dest):
    for pt_journey in pt_journeys.journeys:

        crowflies = [_build_crowfly(pt_journey, **orig), _build_crowfly(pt_journey, **dest)]

        # Filter out empty crowflies
        crowflies = [c for c in crowflies if c]

        pt_journey.sections.extend(crowflies)
        pt_journey.sections.sort(SectionSorter())


def _build_crowfly(pt_journey, requested_obj, mode, places_free_access, fallback_durations, fallback_type):
    fallback_logic = _get_fallback_logic(fallback_type)
    pt_obj = fallback_logic.get_pt_boundaries(pt_journey)

    if pt_obj.uri in places_free_access.odt:
        pt_obj.CopyFrom(requested_obj)
        return None

    section_datetime = fallback_logic.get_section_datetime(pt_journey)
    pt_datetime = fallback_logic.get_pt_bound_datetime(pt_journey)
    fallback_duration = fallback_durations[pt_obj.uri].duration

    crowfly_dt = fallback_logic.get_crowfly_datetime(pt_datetime, fallback_duration)

    crowfly_origin, crowfly_destination = fallback_logic.route_params(requested_obj, pt_obj)
    begin, end = fallback_logic.route_params(crowfly_dt, section_datetime)

    return _create_crowfly(pt_journey, crowfly_origin, crowfly_destination, begin, end, mode)


def _build_fallback(
    requested_obj,
    pt_journeys,
    mode,
    streetnetwork_path_pool,
    obj_accessible_by_crowfly,
    fallback_durations_pool,
    request,
    fallback_type,
):
    accessibles_by_crowfly = obj_accessible_by_crowfly.wait_and_get()
    fallback_durations = fallback_durations_pool.wait_and_get(mode)
    fallback_logic = _get_fallback_logic(fallback_type)

    for pt_journey in pt_journeys.journeys:
        pt_obj = fallback_logic.get_pt_boundaries(pt_journey)

        if requested_obj.uri != pt_obj.uri:

            if pt_obj.uri in accessibles_by_crowfly.odt:
                pt_obj.CopyFrom(requested_obj)
            else:
                # extend the journey with the fallback routing path
                represent_start = fallback_logic.represent_start()
                pt_datetime = fallback_logic.get_pt_bound_datetime(pt_journey)
                fallback_period_extremity = PeriodExtremity(pt_datetime, represent_start)
                orig, dest = fallback_logic.route_params(requested_obj, pt_obj)

                fallback_dp = streetnetwork_path_pool.wait_and_get(
                    orig, dest, mode, fallback_period_extremity, fallback_type, request=request
                )

                if not _is_crowfly_needed(
                    pt_obj.uri, fallback_durations, accessibles_by_crowfly.crowfly, fallback_dp
                ):
                    _replace_crowfly_with_streetnetwork(pt_journey, fallback_dp, fallback_period_extremity)

        pt_journey.sections.sort(SectionSorter())
        fallback_logic.set_fallback_datetime(pt_journey)

    return pt_journeys


def clean_pt_journey_error_or_raise(pt_journeys, has_valid_direct_path_no_pt):
    if pt_journeys and pt_journeys.HasField(b"error"):
        if (
            pt_journeys.error.id == response_pb2.Error.error_id.Value('no_solution')
            and has_valid_direct_path_no_pt
        ):
            pt_journeys.ClearField(b"error")
        else:
            raise PtException(pt_journeys)


def get_max_fallback_duration(request, mode, dp_future):
    """
    By knowing the duration of direct path, we can limit the max duration for proximities by crowfly and fallback
    durations
    :param request:
    :param mode:
    :param dp_future:
    :return:  max_fallback_duration
    """
    # 30 minutes by default
    max_duration = request.get('max_{}_duration_to_pt'.format(mode), 1800)
    dp = dp_future.wait_and_get() if dp_future else None
    dp_duration = dp.journeys[0].durations.total if getattr(dp, 'journeys', None) else max_duration
    return min(max_duration, dp_duration)


def compute_fallback(
    from_obj,
    to_obj,
    streetnetwork_path_pool,
    orig_places_free_access,
    dest_places_free_access,
    request,
    pt_journeys,
):
    """
    Launching fallback computation asynchronously once the pt_journey is finished
    """
    logger = logging.getLogger(__name__)

    for (dep_mode, arr_mode, pt_journey) in pt_journeys:
        logger.debug("completing pt journey that starts with %s and ends with %s", dep_mode, arr_mode)

        if not getattr(pt_journey, "journeys", None):
            continue

        places_free_access = orig_places_free_access.wait_and_get()
        orig_all_free_access = (
            places_free_access.odt | places_free_access.crowfly | places_free_access.free_radius
        )
        places_free_access = dest_places_free_access.wait_and_get()
        dest_all_free_access = (
            places_free_access.odt | places_free_access.crowfly | places_free_access.free_radius
        )

        for journey in pt_journey.journeys:
            # from
            direct_path_type = StreetNetworkPathType.BEGINNING_FALLBACK
            fallback = _get_fallback_logic(direct_path_type)
            pt_orig = fallback.get_pt_boundaries(journey)
            fallback_extremity_dep = PeriodExtremity(journey.departure_date_time, False)
            if from_obj.uri != pt_orig.uri and pt_orig.uri not in orig_all_free_access:
                streetnetwork_path_pool.add_async_request(
                    from_obj, pt_orig, dep_mode, fallback_extremity_dep, request, direct_path_type
                )

            # to
            direct_path_type = StreetNetworkPathType.ENDING_FALLBACK
            fallback = _get_fallback_logic(direct_path_type)
            pt_dest = fallback.get_pt_boundaries(journey)
            fallback_extremity_arr = PeriodExtremity(journey.arrival_date_time, True)
            if to_obj.uri != pt_dest.uri and pt_dest.uri not in dest_all_free_access:
                streetnetwork_path_pool.add_async_request(
                    pt_dest, to_obj, arr_mode, fallback_extremity_arr, request, direct_path_type
                )


def complete_pt_journey(
    requested_orig_obj,
    requested_dest_obj,
    dep_mode,
    arr_mode,
    pt_journeys_,
    streetnetwork_path_pool,
    orig_places_free_access,
    dest_places_free_access,
    orig_fallback_durations_pool,
    dest_fallback_durations_pool,
    request,
):
    """
    We complete the pt journey by adding the beginning fallback and the ending fallback

    :return: pt_journeys: a deeply-copied list of pt_journeys with beginning and ending fallbacks
    """
    logger = logging.getLogger(__name__)

    pt_journeys = copy.deepcopy(pt_journeys_)
    if not getattr(pt_journeys, "journeys", None):
        return pt_journeys

    logger.debug("building pt journey starts with %s and ends with %s", dep_mode, arr_mode)

    pt_journeys = _build_fallback(
        requested_orig_obj,
        pt_journeys,
        dep_mode,
        streetnetwork_path_pool,
        orig_places_free_access,
        orig_fallback_durations_pool,
        request=request,
        fallback_type=StreetNetworkPathType.BEGINNING_FALLBACK,
    )

    pt_journeys = _build_fallback(
        requested_dest_obj,
        pt_journeys,
        arr_mode,
        streetnetwork_path_pool,
        dest_places_free_access,
        dest_fallback_durations_pool,
        request=request,
        fallback_type=StreetNetworkPathType.ENDING_FALLBACK,
    )

    logger.debug("finish building pt journey starts with %s and ends with %s", dep_mode, arr_mode)

    return pt_journeys


def get_entry_point_or_raise(entry_point_obj_future, requested_uri):
    entry_point_obj = entry_point_obj_future.wait_and_get()
    if not entry_point_obj:
        raise EntryPointException(
            "The entry point: {} is not valid".format(requested_uri), response_pb2.Error.unknown_object
        )
    return entry_point_obj


def check_final_results_or_raise(final_results, orig_fallback_durations_pool, dest_fallback_durations_pool):
    if final_results:
        return
    orig_fallback_durations_is_empty = orig_fallback_durations_pool.is_empty()
    dest_fallback_durations_is_empty = dest_fallback_durations_pool.is_empty()
    if orig_fallback_durations_is_empty and dest_fallback_durations_is_empty:
        raise EntryPointException(
            "no origin point nor destination point", response_pb2.Error.no_origin_nor_destination
        )
    if orig_fallback_durations_is_empty:
        raise EntryPointException("no origin point", response_pb2.Error.no_origin)
    if dest_fallback_durations_is_empty:
        raise EntryPointException("no destination point", response_pb2.Error.no_destination)


N_DEG_TO_RAD = 0.01745329238
EARTH_RADIUS_IN_METERS = 6372797.560856


def crowfly_distance_between(start_coord, end_coord):
    lon_arc = (start_coord.lon - end_coord.lon) * N_DEG_TO_RAD
    lon_h = math.sin(lon_arc * 0.5)
    lon_h *= lon_h
    lat_arc = (start_coord.lat - end_coord.lat) * N_DEG_TO_RAD
    lat_h = math.sin(lat_arc * 0.5)
    lat_h *= lat_h
    tmp = math.cos(start_coord.lat * N_DEG_TO_RAD) * math.cos(end_coord.lat * N_DEG_TO_RAD)
    return EARTH_RADIUS_IN_METERS * 2.0 * math.asin(math.sqrt(lat_h + tmp * lon_h))
