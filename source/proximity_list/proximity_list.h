/* Copyright © 2001-2014, Canal TP and/or its affiliates. All rights reserved.
  
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

#pragma once

#include "type/type.h"
#include <vector>
#include <cmath>
#include "nanoflann.h"

namespace navitia { namespace proximitylist {

using type::GeographicalCoord;
struct NotFound : public std::exception{
    NotFound() = default;
    NotFound(const NotFound&) = default;
    NotFound& operator=(const NotFound&) = default;
    virtual ~NotFound() noexcept;
};

/** Définit un indexe spatial qui permet de retrouver les n éléments les plus proches
 *
 * Le template T est le type que l'on souhaite indexer (typiquement un Idx). L'élément sera copié.
 * On rajoute des élements itérativements et on appelle build pour construire l'indexe.
 * L'implémentation est un bête tableau trié par X.
 * On cherche les bornes inf/sup selon X, puis on itère sur les données et on garde les bonnes
 */

template<class T>
struct ProximityList
{
    /// Élement que l'on garde dans le vector 
    struct Item {
        GeographicalCoord coord;
        T element;
        Item(){}
        Item(GeographicalCoord coord, T element) : coord(coord), element(element) {}
        template<class Archive> void serialize(Archive & ar, const unsigned int) {
            ar & coord & element;
        }
    };
    struct Cloud {
        std::vector<Item> items;
        Cloud(const std::vector<Item>& items):items(items){}
        inline size_t kdtree_get_point_count() const { return items.size(); }
        inline double kdtree_distance(const double *p1, const size_t idx_p2, size_t /*size*/) const
        {
            GeographicalCoord tmp_coord{p1[0], p1[1]};
            static const double DEG_TO_RAD = 0.0174532925199432958;
            const double rad = tmp_coord.lat() * DEG_TO_RAD;

            /*const double coslat = ::cos(tmp_coord.lat() * DEG_TO_RAD)*/;
            const double coslat  = 1 - rad*rad/2. +  rad*rad* rad*rad/24.- rad*rad* rad*rad*rad*rad/720.;
            return tmp_coord.approx_sqr_distance(items[idx_p2].coord, coslat);
        }
        inline double kdtree_get_pt(const size_t idx, int dim) const
        {
            if (dim==0) return items[idx].coord.lon();
            else return items[idx].coord.lat();
        }
        template <class BBOX>
        bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }
    };
    typedef nanoflann::KDTreeSingleIndexAdaptor<
            nanoflann::L2_Simple_Adaptor<double, Cloud> ,
            Cloud,
            2 /* dim */
        > my_kd_tree_t;
    mutable std::unique_ptr<my_kd_tree_t> my_tree=nullptr;
    mutable std::unique_ptr<Cloud> cloud = nullptr;
    /// Contient toutes les coordonnées de manière à trouver rapidement
    std::vector<Item> items;

    /// Rajoute un nouvel élément. Attention, il faut appeler build avant de pouvoir utiliser la structure
    void add(GeographicalCoord coord, T element){
        items.push_back(Item(coord,element));
    }
    void clear(){
        items.clear();
    }

    /// Construit l'indexe
    void build();

    /// Retourne tous les éléments dans un rayon de x mètres
    std::vector< std::pair<T, GeographicalCoord> > find_within(GeographicalCoord coord, double distance = 500)  const ;

    /// Fonction de confort pour retrouver l'élément le plus proche dans l'indexe
    T find_nearest(double lon, double lat) const {
        return find_nearest(GeographicalCoord(lon, lat));
    }

    /// Retourne l'élément le plus proche dans tout l'indexe
    T find_nearest(GeographicalCoord coord, double max_dist = 500) const {
        auto temp = find_within(coord, max_dist);
        if(temp.empty())
            throw NotFound();
        else

            return temp.front().first;
    }

    /** Fonction qui permet de sérialiser (aka binariser la structure de données
      *
      * Elle est appelée par boost et pas directement
      */
    template<class Archive> void serialize(Archive & ar, const unsigned int) {
        ar & items;
    }

};

}} // namespace navitia::proximitylist
