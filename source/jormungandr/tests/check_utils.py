# Copyright (c) 2001-2014, Canal TP and/or its affiliates. All rights reserved.
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

from collections import deque, defaultdict
from nose.tools import *
import json
from navitiacommon import request_pb2, response_pb2
from datetime import datetime
import logging
import re
from shapely.geometry import shape
import sys


"""
Some small functions to check the service responses
"""


def check_url(tester, url, might_have_additional_args=False, **kwargs):
    """
    Test url status code to 200 and if valid format response as json
    if might_have_additional_args is set to True,
        we just don't want a 404 error (the url might be in error because of mandatory params not provided)
    else
        we don't want an error on the url
    """
    response = tester.get(url, **kwargs)

    assert response, "response for url {} is null".format(url)
    if might_have_additional_args:
        assert response.status_code != 404, "unreachable url {}"\
            .format(json.dumps(json.loads(response.data), indent=2))
    else:
        eq_(response.status_code, 200, "invalid return code, response : {}"
            .format(json.dumps(json.loads(response.data, encoding='utf-8'), indent=2)))
    return response


def get_not_null(dict, field):
    assert field in dict

    val = dict[field]
    if type(val) == bool:
        return val  # no check for booleans
    if type(val) == int:
        return val  # no check for integer

    assert val, "value of field {} is null".format(field)
    return val


days_regexp = re.compile("^(0|1){366}$")


def is_valid_days(days):
    m = days_regexp.match(days)
    return m is not None


version_number_regexp = re.compile("v[0-9]+\.[0-9]+\.[0-9]+[-.*]?")


def is_valid_navitia_version_number(str):
    """
    check that the version number is valid
    it must contains at least v{major}.{minor}.{hotfix}
    it can also contains the git sha1 at the end
    >>> is_valid_navitia_version_number("v1.12.126")
    True
    >>> is_valid_navitia_version_number("v1.3.1-73-g4c7524b")
    True
    >>> is_valid_navitia_version_number("1.12.126")
    Traceback (most recent call last):
    AssertionError
    >>> is_valid_navitia_version_number("v12.126-73-g4c7524b")
    Traceback (most recent call last):
    AssertionError
    """
    m = version_number_regexp.match(str)
    assert m
    return True


def get_valid_datetime(str, possible_errors=False):
    """
    Check is the string is a valid date and return it
    if possible_errors, the string might be equals to
    "not-a-date-time"

    >>> get_valid_datetime("bob")
    Traceback (most recent call last):
    AssertionError
    >>> get_valid_datetime("")
    Traceback (most recent call last):
    AssertionError
    >>> get_valid_datetime("20123101T215030")  # month is badly set
    Traceback (most recent call last):
    AssertionError
    >>> get_valid_datetime("20123101T215030", possible_errors=True)
    Traceback (most recent call last):
    AssertionError
    >>> get_valid_datetime("not-a-date-time", possible_errors=True)

    >>> get_valid_datetime("not-a-date-time")
    Traceback (most recent call last):
    AssertionError
    >>> get_valid_datetime("20120131T215030")
    datetime.datetime(2012, 1, 31, 21, 50, 30)
    """
    assert str

    try:
        return datetime.strptime(str, "%Y%m%dT%H%M%S")
    except ValueError:
        if possible_errors:
            assert str == "not-a-date-time"
            return None
        logging.error("string '{}' is no valid datetime".format(str))
        assert False


def get_valid_time(str):
    """
    Check is the string is a valid time and return it
    >>> get_valid_time("bob")
    Traceback (most recent call last):
    AssertionError
    >>> get_valid_time("")
    Traceback (most recent call last):
    AssertionError
    >>> get_valid_time("20120131T215030")  # it's a datetime, not valid
    Traceback (most recent call last):
    AssertionError
    >>> get_valid_time("215030")  #time is HHMMSS
    datetime.datetime(1900, 1, 1, 21, 50, 30)
    >>> get_valid_time("501230")  # MMHHSS, not valid
    Traceback (most recent call last):
    AssertionError
    """
    assert str

    try:
        #AD:we use a datetime anyway because I don't know what to use instead
        return datetime.strptime(str, "%H%M%S")
    except ValueError:
        logging.error("string '{}' is no valid time".format(str))
        assert False


def is_valid_date(str):
    """
    Check is the string is a valid date
    >>> is_valid_date("bob")
    False
    >>> is_valid_date("")
    Traceback (most recent call last):
    AssertionError
    >>> is_valid_date("20123101")  # month is badly set
    False
    >>> is_valid_date("20120131")
    True
    """
    assert str

    try:
        datetime.strptime(str, "%Y%m%d")
    except ValueError:
        logging.error("string '{}' is no valid date".format(str))
        return False
    return True


def is_valid_bool(str):
    if type(str) is bool:
        return True

    assert str
    #else check as string
    lower = str.lower()
    return lower == "true" or lower == "false"


def get_valid_float(str):
    if type(str) is float:
        return str

    try:
        return float(str)
    except ValueError:
        assert "cannot convert {} to float".format(str)


def get_valid_int(str):
    assert str != ""
    if type(str) is int:
        return str

    try:
        return int(str)
    except ValueError:
        assert False


def is_valid_lat(str):
    lat = get_valid_float(str)

    assert -90.0 <= lat <= 90.0, "lat should be between -90 and 90"


def is_valid_lon(str):
    lat = get_valid_float(str)

    assert 180.0 >= lat >= -180.0, "lon should be between -180 and 180"


def is_valid_coord(coord):
    lat = get_not_null(coord, "lat")
    lon = get_not_null(coord, "lon")
    is_valid_lat(lat)
    is_valid_lon(lon)


def get_links_dict(response):
    """
    get links as dict ordered by 'rel' or 'type"
    """
    raw_links = get_not_null(response, "links")

    #create a dict with the 'rel' field as key
    links = {link.get('rel', link.get('type', None)): link for link in raw_links}

    return links


def check_links(object, tester, href_mandatory=True):
    """
    get the links as dict ordered by 'rel' and check:
     - all links must have the attributes:
       * 'internal' --> optional but must be a boolean
       * 'href' --> valid url if not templated, empty if internal
       * 'rel' --> not empty if internal
       * 'title' --> optional
       * 'templated' --> optional but must be a boolean
       * 'type' --> not empty
    """
    links = get_links_dict(object)

    for link_name, link in links.iteritems():
        def get_bool(name):
            """ give boolean if in dict, else False"""
            if name in link:
                assert is_valid_bool(link[name])
                if bool(link[name]):
                    return True
            return False
        internal = get_bool('internal')
        templated = get_bool('templated')

        if href_mandatory:
            if not internal:
                assert 'href' in link, "no href in link"

            if not templated and not internal:
                #we check that the url is valid
                assert check_url(tester, link['href'].replace('http://localhost', ''),
                                 might_have_additional_args=False), "href's link must be a valid url"

        if internal:
            assert 'rel' in link
            assert link['rel']

        assert 'type' in link
        assert link['type']

    return links


def check_internal_links(response, tester):
    """
    We want to check that all 'internal' link are correctly linked to an element in the response

    for that we first collect all internal link
    then iterate on all node and remove a link if we find a node with
     * a name equals to link.'rel'
     * an id equals to link.'id'

     At the end the internal link list must be empty
    """
    from jormungandr import utils  #import late not to load it before updating the conf for the tests

    internal_links_id = set()
    internal_link_types = set()  # set with the types we look for

    def add_link_visitor(name, val):
        if val and name == 'links':
            if 'internal' in val and bool(val['internal']):
                internal_links_id.add(val['id'])
                internal_link_types.add(val['rel'])

    utils.walk_dict(response, add_link_visitor)

    def check_node(name, val):

        if name in internal_link_types:

            if 'id' in val and val['id'] in internal_links_id:
                #found one ref, we can remove the link
                internal_links_id.remove(val['id'])

    utils.walk_dict(response, check_node)

    assert not internal_links_id, "cannot find correct ref for internal links : {}".\
        format([lid for lid in internal_links_id])


class unique_dict(dict):
    """
    We often have to check that a set of values are uniq, this container is there to do the job

    >>> d = unique_dict('id')
    >>> d['bob'] = 1
    >>> d['bobette'] = 1
    >>> d['bob'] = 2
    Traceback (most recent call last):
        ...
    AssertionError: the id if must be unique, but 'bob' is not
    """

    def __init__(self, key_name):
        self.key_name = key_name

    def __setitem__(self, key, value):
        assert not key in self, \
            "the {} if must be unique, but '{}' is not".format(self.key_name, key)
        dict.__setitem__(self, key, value)


def query_from_str(str):
    """
    for convenience, convert a url to a dict

    >>> query_from_str("toto/tata?bob=toto&bobette=tata&bobinos=tutu")
    {'bobette': 'tata', 'bobinos': 'tutu', 'bob': 'toto'}
    >>> query_from_str("toto/tata?bob=toto&bob=tata&bob=titi&bob=tata&bobinos=tutu")
    {'bobinos': 'tutu', 'bob': ['toto', 'tata', 'titi', 'tata']}

    Note: the query can be encoded, so the split it either on the encoded or the decoded value
    """
    query = {}
    last_elt = str.split('?' if '?' in str else '%3F')[-1]

    for s in last_elt.split('&' if '&' in last_elt else '%26'):
        k, v = s.split("=" if '=' in s else '%3D')

        if k in query:
            old_val = query[k]
            if isinstance(old_val, list):
                old_val.append(v)
            else:
                query[k] = [old_val, v]
        else:
            query[k] = v

    return query


def is_valid_feed_publisher(feed_publisher):
    get_not_null(feed_publisher, 'id')
    get_not_null(feed_publisher, 'name')
    get_not_null(feed_publisher, 'license')
    get_not_null(feed_publisher, 'url')


def is_valid_journey_response(response, tester, query_str):
    if isinstance(query_str, basestring):
        query_dict = query_from_str(query_str)
    else:
        query_dict = query_str

    journeys = get_not_null(response, "journeys")

    all_sections = unique_dict('id')
    assert len(journeys) > 0, "we must at least have one journey"
    for j in journeys:
        is_valid_journey(j, tester, query_dict)

        for s in j['sections']:
            all_sections[s['id']] = s

    # check the fare section
    # the fares must be structurally valid and all link to sections must be ok
    all_tickets = unique_dict('id')
    fares = response['tickets']
    for f in fares:
        is_valid_ticket(f, tester)
        all_tickets[f['id']] = f

    check_internal_links(response, tester)

    #check other links
    check_links(response, tester)

    # more checks on links, we want the prev/next/first/last,
    # to have forwarded all params, (and the time must be right)
    journeys_links = get_links_dict(response)

    for l in ["prev", "next", "first", "last"]:
        assert l in journeys_links
        url = journeys_links[l]['href']

        additional_args = query_from_str(url)
        for k, v in additional_args.iteritems():
            if k == 'datetime':
                #TODO check datetime
                continue
            if k == 'datetime_represents':
                query_dt_rep = query_dict.get('datetime_represents', 'departure')
                if l in ['prev', 'last']:
                    #the datetime_represents is negated
                    if query_dt_rep == 'departure':
                        assert v == 'arrival'
                    else:
                        assert v == 'departure'
                else:
                    query_dt_rep == v

                continue

            eq_(query_dict[k], v)

    feed_publishers = get_not_null(response, "feed_publishers")
    for feed_publisher in feed_publishers:
        is_valid_feed_publisher(feed_publisher)


def is_valid_journey(journey, tester, query):
    arrival = get_valid_datetime(journey['arrival_date_time'])
    departure = get_valid_datetime(journey['departure_date_time'])
    request = get_valid_datetime(journey['requested_date_time'])

    assert arrival >= departure
    # test if duration time is consistent with arrival and departure
    # as we sometimes loose a second in rounding section duration, tolerance is added
    assert (arrival - departure).seconds - journey['duration'] <= len(journey['sections']) - 1

    if 'datetime_represents' not in query or query['datetime_represents'] == "departure":
        #for 'departure after' query, the departure must be... after \o/
        assert departure >= request
    else:
        assert arrival <= request

    #we want to test that all departure match de previous section arrival
    last_arrival = departure
    for s in journey['sections']:
        is_valid_section(s, query)
        section_departure = get_valid_datetime(s['departure_date_time'])
        assert (section_departure - last_arrival).seconds <= 1  # there cannot be more than one second between the 2

        last_arrival = get_valid_datetime(s['arrival_date_time'])

        # test if geojson is valid
        g = s.get('geojson')
        g is None or shape(g)

    assert last_arrival == arrival
    assert get_valid_datetime(journey['sections'][-1]['arrival_date_time']) == last_arrival


def is_valid_section(section, query):
    arrival = get_valid_datetime(section['arrival_date_time'])
    departure = get_valid_datetime(section['departure_date_time'])

    assert (arrival - departure).seconds == section['duration']

    assert section['type']  # type cannot be empty

    #for street network section, we must have a valid path
    if section['type'] == 'street_network':
        assert section['mode']  # mode cannot be empty for street network
        total_duration = 0
        for p in section['path']:
            assert get_valid_int(p['length']) >= 0
            assert -180 <= get_valid_int(p['direction']) <= 180  # direction is an angle
            #No constraint on name, it can be empty
            dur = get_valid_int(p['duration'])
            assert dur >= 0
            total_duration += dur

        assert total_duration == section['duration']

    #TODO check geojson
    #TODO check stop_date_times
    #TODO check from/to


def is_valid_ticket(ticket, tester):
    found = get_not_null(ticket, 'found')
    assert is_valid_bool(found)

    get_not_null(ticket, 'id')
    get_not_null(ticket, 'name')
    cost = get_not_null(ticket, 'cost')
    if found:
        #for found ticket, we must have a non empty currency
        get_not_null(cost, 'currency')

    get_valid_float(get_not_null(cost, 'value'))

    check_links(ticket, tester)


def is_valid_stop_area(stop_area, depth_check=1):
    """
    check the structure of a stop area
    """
    get_not_null(stop_area, "name")
    coord = get_not_null(stop_area, "coord")
    is_valid_label(get_not_null(stop_area, "label"))
    is_valid_coord(coord)

    for c in stop_area.get('comments', []):
        is_valid_comment(c)


def is_valid_stop_point(stop_point, depth_check=1):
    """
    check the structure of a stop point
    """

    get_not_null(stop_point, "name")
    is_valid_label(get_not_null(stop_point, "label"))
    coord = get_not_null(stop_point, "coord")
    is_valid_coord(coord)

    for c in stop_point.get('comments', []):
        is_valid_comment(c)

    if depth_check > 0:
        is_valid_stop_area(get_not_null(stop_point, "stop_area"), depth_check-1)
    else:
        assert "stop_area" not in stop_point


def is_valid_route(route, depth_check=1):
    get_not_null(route, "name")
    is_valid_bool(get_not_null(route, "is_frequence"))

    direction = get_not_null(route, "direction")
    is_valid_place(direction, depth_check - 1)
    #the direction of the route must always be a stop point
    assert get_not_null(direction, "embedded_type") == "stop_area"
    is_valid_stop_area(get_not_null(direction, "stop_area"), depth_check - 1)

    if depth_check > 0:
        is_valid_line(get_not_null(route, "line"), depth_check - 1)
    else:
        assert 'line' not in route

    for c in route.get('comments', []):
        is_valid_comment(c)

    # test if geojson is valid
    g = route.get('geojson')
    g is None or shape(g) #TODO check length


def is_valid_company(company, depth_check=1):
    get_not_null(company, "name")
    get_not_null(company, "id")


def is_valid_physical_mode(physical_mode, depth_check=1):
    get_not_null(physical_mode, "name")
    get_not_null(physical_mode, "id")


def is_valid_line(line, depth_check=1):
    get_not_null(line, "name")
    get_not_null(line, "id")

    for c in line.get('comments', []):
        is_valid_comment(c)

    if depth_check > 0:
        is_valid_network(get_not_null(line, 'network'), depth_check - 1)

        routes = get_not_null(line, 'routes')
        for r in routes:
            is_valid_route(r, depth_check - 1)
    else:
        assert 'network' not in line
        assert 'routes' not in line

    # test if geojson is valid
    g = line.get('geojson')
    g is None or shape(g) #TODO check length

def is_valid_line_group(line_group, depth_check=1):
    get_not_null(line_group, "name")
    get_not_null(line_group, "id")

    if depth_check > 0:
        # the main_line is always displayed with a depth of 0 to reduce duplicated informations
        is_valid_line(get_not_null(line_group, "main_line"), 0)
        for l in line_group.get('lines', []):
            is_valid_line(l, depth_check - 1)


def is_valid_poi(poi, depth_check=1):
    get_not_null(poi, 'name')
    poi_type = get_not_null(poi, 'poi_type')
    get_not_null(poi_type, 'id')
    get_not_null(poi_type, 'name')
    get_not_null(poi, 'label')
    get_not_null(poi, 'id')
    is_valid_coord(get_not_null(poi, 'coord'))
    for admin in get_not_null(poi, 'administrative_regions'):
        is_valid_admin(admin, depth_check-1)
    is_valid_address(get_not_null(poi, 'address'), depth_check - 1)


def is_valid_admin(admin, depth_check=1):
    if depth_check < 0:
        return
    get_not_null(admin, 'insee')
    name = get_not_null(admin, 'name')
    zip_code = get_not_null(admin, 'zip_code')
    lbl = get_not_null(admin, 'label')
    is_valid_label(lbl)
    assert name in lbl and zip_code in lbl  # so name of the admin and it's zip code must be in the label

    get_not_null(admin, 'id')
    get_valid_int(get_not_null(admin, 'level'))
    is_valid_coord(get_not_null(admin, 'coord'))


def is_valid_codes(codes):
    for code in codes:
        get_not_null(code, "type")
        get_not_null(code, "value")


def is_valid_places(places, depth_check=1):
    for p in places:
        is_valid_place(p, depth_check)


def is_valid_place(place, depth_check=1):
    if depth_check < 0:
        return
    n = get_not_null(place, "name")
    get_not_null(place, "id")
    type = get_not_null(place, "embedded_type")
    if type == "address":
        address = get_not_null(place, "address")
        is_valid_address(address, depth_check)
    elif type == "stop_area":
        stop_area = get_not_null(place, "stop_area")
        is_valid_stop_area(stop_area, depth_check)

        #for stops name should be the label
        is_valid_label(n)
        assert stop_area['label'] == n
    elif type == "stop_point":
        stop_point = get_not_null(place, "stop_point")
        is_valid_stop_point(stop_point, depth_check)

        is_valid_label(n)
        assert stop_point['label'] == n
    elif type == "poi":
        is_valid_poi(get_not_null(place, "poi"), depth_check)
    else:
        assert(False, "invalid type")


def is_valid_pt_objects_response(response, depth_check=1):
    for pt_obj in get_not_null(response, 'pt_objects'):
        is_valid_pt_object(pt_obj, depth_check)


def is_valid_pt_object(pt_object, depth_check=1):
    n = get_not_null(pt_object, "name")
    get_not_null(pt_object, "id")
    get_not_null(pt_object, "quality")
    pt_obj_type = get_not_null(pt_object, "embedded_type")

    assert pt_obj_type in ('line',
                           'route',
                           'network',
                           'commercial_mode',
                           'admin',
                           'vehicle_journey',
                           'calendar',
                           'company',
                           "stop_area",
                           "stop_point",
                           "poi",
                           "address")

    # if it's a line, it should pass the 'is_valid_line' test,
    # if it's a stop_area, it should pass the 'is_valid_stop_area' test ...
    check_method_to_call = getattr(sys.modules[__name__], 'is_valid_' + pt_obj_type)
    check_method_to_call(get_not_null(pt_object, pt_obj_type), depth_check)

    # check the pt_object label
    if pt_obj_type in ('stop_area', 'stop_point'):
        #for stops name should be the label
        assert get_not_null(pt_object, pt_obj_type)['label'] == n
    if pt_obj_type == 'line':
        # the line network, commercial_mode, code and name should be in the label
        check_embedded_line_label(n, pt_obj_type['line'], depth_check)
    if pt_obj_type == 'route':
        # the line network, commercial_mode, code and name should be in the label
        check_embedded_route_label(n, pt_obj_type['route'], depth_check)


def check_embedded_line_label(label, line, depth_check):
    is_valid_label(label)
    #the label must contains
    # the network name, the commercial mode name, the line id, and the line name
    if depth_check > 0:
        assert get_not_null(line, 'commercial_mode')['name'] in label
        assert get_not_null(line, 'network')['name'] in label
    assert get_not_null(line, 'code') in label
    assert get_not_null(line, 'name') in label


def check_embedded_route_label(label, route, depth_check):
    is_valid_label(label)
    #the label must contains
    # the line's network name, the commercial mode name, the line id, and the route name
    assert get_not_null(route, 'name') in label
    if depth_check > 0:
        line = get_not_null(route, 'line')
        assert get_not_null(line, 'code') in label

        if depth_check > 1:
            assert get_not_null(line, 'commercial_mode')['name'] in label
            assert get_not_null(line, 'network')['name'] in label


def is_valid_address(address, depth_check=1):
    if depth_check < 0:
        return
    id = get_not_null(address, "id")
    lon, lat = id.split(';')
    is_valid_lon(lon)
    is_valid_lat(lat)
    get_not_null(address, "house_number")
    get_not_null(address, "name")
    if depth_check >= 1:
        for admin in get_not_null(address, "administrative_regions"):
            is_valid_admin(admin, depth_check)
    coord = get_not_null(address, "coord")
    is_valid_coord(coord)


def is_valid_validity_pattern(validity_pattern, depth_check=1):
    beginning_date = get_not_null(validity_pattern, "beginning_date")
    assert is_valid_date(beginning_date)

    days = get_not_null(validity_pattern, "days")
    assert is_valid_days(days)


def is_valid_network(network, depth_check=1):
    get_not_null(network, "id")
    get_not_null(network, "name")


def is_valid_vehicle_journey(vj, depth_check=1):
    if depth_check < 0:
        return
    get_not_null(vj, "id")
    get_not_null(vj, "name")

    for c in vj.get('comments', []):
        is_valid_comment(c)

    if depth_check > 0:
        is_valid_journey_pattern(get_not_null(vj, 'journey_pattern'), depth_check=depth_check-1)
        is_valid_validity_pattern(get_not_null(vj, 'validity_pattern'), depth_check=depth_check-1)

        stoptimes = get_not_null(vj, 'stop_times')

        for st in stoptimes:
            get_valid_time(get_not_null(st, 'arrival_time'))
            get_valid_time(get_not_null(st, 'departure_time'))

            if depth_check > 1:
                #with depth > 1 (we are already in the stoptime nested object), we don't want jpp
                is_valid_journey_pattern_point(get_not_null(st, 'journey_pattern_point'), depth_check - 2)
            else:
                assert 'journey_pattern_point' not in st
    else:
        #with depth = 0, we don't want the stop times, the jp, vp, ...
        assert 'stop_times' not in vj
        assert 'journey_pattern' not in vj
        assert 'validity_pattern' not in vj


def is_valid_journey_pattern(jp, depth_check=1):
    if depth_check < 0:
        return
    get_not_null(jp, "id")
    get_not_null(jp, "name")


def is_valid_journey_pattern_point(jpp, depth_check=1):
    get_not_null(jpp, "id")
    if depth_check > 0:
        is_valid_stop_point(get_not_null(jpp, 'stop_point'), depth_check=depth_check - 1)
    else:
        assert 'stop_point' not in jpp


def is_valid_comment(comment):
    get_not_null(comment, 'type')
    get_not_null(comment, 'value')


def is_valid_region_status(status):
    get_not_null(status, 'status')
    get_valid_int(get_not_null(status, 'data_version'))
    get_valid_int(get_not_null(status, 'nb_threads'))
    is_valid_bool(get_not_null(status, 'last_load_status'))
    is_valid_bool(get_not_null(status, 'is_connected_to_rabbitmq'))
    is_valid_date(get_not_null(status, 'end_production_date'))
    is_valid_date(get_not_null(status, 'start_production_date'))
    get_valid_datetime(get_not_null(status, 'last_load_at'), possible_errors=True)
    get_valid_datetime(get_not_null(status, 'publication_date'), possible_errors=True)


# for () are mandatory for the label even if is reality it is not
# (for example if the admin has no post code or a stop no admin)
# This make the test data a bit more difficult to create, but that way we can check the label creation
label_regexp = re.compile(".* \(.*\)")


def is_valid_label(label):
    m = label_regexp.match(label)
    return m is not None


def get_disruptions(obj, response):
    """
    unref disruption links are return the list of disruptions
    """
    all_disruptions = {d['id']: d for d in response['disruptions']}

    if 'links' not in obj:
        return None
    return [all_disruptions[d['id']] for d in obj['links'] if d['type'] == 'disruption']


def is_valid_disruption(disruption):
    get_not_null(disruption, 'id')
    get_not_null(disruption, 'disruption_id')
    s = get_not_null(disruption, 'severity')
    get_not_null(s, 'name')
    get_not_null(s, 'color')
    get_not_null(s, 'effect')
    msg = get_not_null(disruption, 'messages')
    assert len(msg) > 0
    for m in msg:
        get_not_null(m, "text")
        channel = get_not_null(m, 'channel')
        get_not_null(channel, "content_type")
        get_not_null(channel, "id")
        get_not_null(channel, "name")

s_coord = "0.0000898312;0.0000898312"  # coordinate of S in the dataset
r_coord = "0.00188646;0.00071865"  # coordinate of R in the dataset

#default journey query used in various test
journey_basic_query = "journeys?from={from_coord}&to={to_coord}&datetime={datetime}"\
    .format(from_coord=s_coord, to_coord=r_coord, datetime="20120614T080000")

def get_all_disruptions(elem, response):
    """
    return a map with the disruption id as key and the list of disruption + impacted object as value for a item of the response
    """
    disruption_by_obj = defaultdict(list)

    all_disruptions = {d['id']: d for d in response['disruptions']}

    def disruptions_filler(_, obj):
        try:
            if 'links' not in obj:
                return
        except TypeError:
            return

        real_disruptions = [all_disruptions[d['id']] for d in obj['links'] if d['type'] == 'disruption']

        for d in real_disruptions:
            disruption_by_obj[d['id']].append((d, obj))

    #we import utils here else it will import jormungandr too early in the test
    from jormungandr import utils
    utils.walk_dict(elem, disruptions_filler)

    return disruption_by_obj


def is_valid_stop_date_time(stop_date_time):
    get_not_null(stop_date_time, 'arrival_date_time')
    assert get_valid_datetime(stop_date_time['arrival_date_time'])
    get_not_null(stop_date_time, 'departure_date_time')
    assert get_valid_datetime(stop_date_time['departure_date_time'])

def get_interpreter():
    """
    Return the current interpreter(python, pypy.... etc.)
    It should only be used in unit tests.
    :return: str
    """
    # sys.executable is the whole path of the interpreter
    exe = sys.executable

    if 'python' in exe:
        return 'python'
    if 'pypy' in exe:
        return 'pypy'
