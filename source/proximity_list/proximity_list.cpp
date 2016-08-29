#include "proximity_list.h"


namespace navitia { namespace proximitylist {

NotFound::~NotFound() noexcept {}

template<typename T>
void ProximityList<T>::build(){
    //std::sort(items.begin(), items.end(), [](const Item & a, const Item & b){return a.coord < b.coord;});
    cloud = std::make_unique<Cloud>(items);
    my_tree = std::make_unique<my_kd_tree_t>(2 /*dim*/, *cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
    my_tree->buildIndex();
}

/// Retourne tous les éléments dans un rayon de x mètres
template<typename T>
std::vector< std::pair<T, GeographicalCoord> > ProximityList<T>::find_within(GeographicalCoord coord, double distance) const {
    /*
    double distance_degree = distance / 111320;

    static const double DEG_TO_RAD = 0.0174532925199432958;

    const double coslat = ::cos(coord.lat() * DEG_TO_RAD);

    const double max_lat = coord.lat() + distance_degree;
    const double min_lat = coord.lat() - distance_degree;

    auto begin = std::lower_bound(items.begin(), items.end(), coord.lon() - distance_degree / coslat, [](const Item & i, double min){return i.coord.lon() < min;});
    auto end = std::upper_bound(begin, items.end(), coord.lon() + distance_degree / coslat, [](double max, const Item & i){return max < i.coord.lon();});

    std::vector< std::pair<T, GeographicalCoord> > result;

    double max_dist = distance * distance;
    for(; begin != end; ++begin){
        if( min_lat < begin->coord.lat() && begin->coord.lat() < max_lat && begin->coord.approx_sqr_distance(coord, coslat) <= max_dist)
            result.push_back(std::make_pair(begin->element, begin->coord));
    }
    std::sort(result.begin(), result.end(), [&coord, &coslat](const std::pair<T, GeographicalCoord> & a, const std::pair<T, GeographicalCoord> & b){return a.second.approx_sqr_distance(coord, coslat) < b.second.approx_sqr_distance(coord, coslat);});
    */
    std::vector<std::pair<size_t, double>> ret_matches;
    const double query_pt[2] = { coord.lon(), coord.lat()};
    nanoflann::SearchParams params;
    size_t nMatches = 0;
    if (!my_tree) {
        cloud = std::make_unique<Cloud>(items);
        my_tree = std::make_unique<my_kd_tree_t>(2 /*dim*/, *cloud, nanoflann::KDTreeSingleIndexAdaptorParams(128 /* max leaf */));
        my_tree->buildIndex();
    }
    if (my_tree) {
         nMatches=my_tree->radiusSearch(query_pt, 500*500, ret_matches, params);
    }
    std::vector< std::pair<T, GeographicalCoord> > result;
    result.reserve(nMatches);
    for (auto& p: ret_matches) {
        auto item = cloud->items[p.first];
        result.push_back({item.element, item.coord});
    }
    return result;
}
template class ProximityList<unsigned int>;
template class ProximityList<unsigned long>;
}} // namespace navitia::proximitylist
