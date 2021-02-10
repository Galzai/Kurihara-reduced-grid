#pragma once

#include "CoordinatesMath.h"
#include "GISNamedTypes.h"
#include <typeinfo>
#include <typeindex>
#include <concepts>
#include <ranges>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <stack>

// Two concepts required by the below functions:

template<class Me, class Other> concept derived_or_same =
    std::same_as<Me, Other> || std::derived_from<Me, Other>;

template<class Me, class Other> concept concrete_derived_or_same =
	!std::is_abstract_v<Me> && derived_or_same<Me, Other>;

// class Grid
template<typename Entity, std::size_t num_rows> requires (num_rows > 0)
class Grid { 
    struct GridRow;
    class Cell;
    friend class FooTest;

    Latitude m_cellHeightDeg;
    Meters m_cellHeightMeters;
    std::size_t m_numCells;
    mutable std::vector<GridRow> m_gridCells;
    std::list<Cell> m_cells;

    // We need an epsilon for float comparison
    static constexpr double Epsilon = 0.00002;

    // Simple struct for compact representation of row
    struct GridRow{
        Longitude cellWidthDeg;
        std::vector<Cell*> rowCells;
        GridRow(Longitude cellWidthVal, std::vector<Cell*> rowCellsVal) : cellWidthDeg(cellWidthVal), rowCells(rowCellsVal){};
    };
    
    // private inner class Cell
    class Cell {
        mutable std::unordered_map<std::type_index, std::vector<Entity*>> m_entitiesByType;
        mutable std::vector<Entity*> m_allEntities;

    public:
        // Add an entity of concrete type to cell by type and to all
        template<typename ActualT> requires derived_or_same<ActualT, Entity>
        void add(ActualT& e) {
            m_entitiesByType[typeid(ActualT)].push_back(&e);
            m_allEntities.push_back(&e);
        }  
        // Cell’s Constructors and Assignment:
        // There are no constructors that are required by the API. You should decide whether to implement, block or rely on the default for copy and move constructors, as well as for the assignment operators.
        Cell(){};
        // Cell’s Destructor:
        // You should decide whether there is a need to implement a user defined destructor.

        // (Cell::A1) Getting Entities from a Cell:
        // returns all entities in that cell, which return true for the PredicateFunc.
        template<typename PredicateFunc>
        std::vector<Entity*> getEntities(PredicateFunc&& func) const {
            std::vector<Entity*> filteredEntities;
            std::copy_if(m_allEntities.begin(), m_allEntities.end(), std::back_inserter(filteredEntities), 
            [func](Entity *entity){ return func(*entity);} );
            return filteredEntities;
        }

        // (Cell::B1) 
        // Getting specific type of Entities from a Cell:
        // returns all entities of type ActualT which return true for the PredicateFunc.
        template<typename ActualT, typename PredicateFunc> 
        requires concrete_derived_or_same<ActualT, Entity>
        std::vector<ActualT*> getEntities(PredicateFunc&& func) const {
            std::vector<ActualT*> filteredEntities;
            if(m_entitiesByType.find(typeid(ActualT)) == m_entitiesByType.end()) return filteredEntities;
            // copy untill no more to copy
            for(auto p_entity : m_entitiesByType.at(typeid(ActualT)))
            {
                ActualT *p_actualEntity = dynamic_cast<ActualT*>(p_entity);
                if(func(*p_actualEntity)){
                    filteredEntities.push_back(p_actualEntity);
                }
            }
            return filteredEntities;
        }

        // (Cell::B2) 
        // same as (Cell::B1) above but with a limit on the number of returned entities (up to limit entities).
        template<typename ActualT, typename PredicateFunc> requires concrete_derived_or_same<ActualT, Entity>
        std::vector<ActualT*> getEntities(PredicateFunc&& func, std::size_t limit) const {
            std::vector<ActualT*> filteredEntities;
            std::size_t curCopied = 0;
            if(m_entitiesByType.find(typeid(ActualT)) == m_entitiesByType.end()) return filteredEntities;
            // copy untill full or untill no more to copy
            for(auto p_entity : m_entitiesByType.at(typeid(ActualT)))
            {
                ActualT *p_actualEntity = dynamic_cast<ActualT*>(p_entity);
                if(func(*p_actualEntity)){
                    filteredEntities.push_back(p_actualEntity);
                    if(++curCopied == limit) break;
                }
            }
            return filteredEntities;
        }

        // (Cell::B3) 
        // returns a range of all entities of type ActualT.
        // Complexity of this function is O(1)
        // This function returns a view that is updated “behind the scene” automatically in case additional objects of type ActualT are added to this Cell via the Grid. Order of entities in the view shall be the same as the order of their insertion to the grid.
        // Iterating over the returned view generates pointers with the correct type, i.e. ActualT*.
        template<typename ActualT> requires concrete_derived_or_same<ActualT, Entity>
        std::ranges::sized_range auto getEntitiesView() const {
            //  Helper class as given by instructor, used to to iterate by type in O(1)
            class ViewToVectorOfPointers {
                std::vector<Entity*>& vec;
            public:
            // Struct for the iterator
                struct PointersViewIterator: public std::vector<Entity*>::iterator {

                    // Constructors
                    PointersViewIterator() {}
                    PointersViewIterator(std::vector<Entity*>::iterator itr): std::vector<Entity*>::iterator(itr) {}

                    // Required operators
                    ActualT* operator*() { return static_cast<ActualT*>(this->std::vector<Entity*>::iterator::operator*()); }
                    const ActualT* operator*() const { return static_cast<const ActualT*>(this->std::vector<Entity*>::iterator::operator*()); }
                    PointersViewIterator& operator++() {
                        this->std::vector<Entity*>::iterator::operator++();
                        return *this;
                    }
                    PointersViewIterator operator++(int) {
                        auto old = *this;
                        ++(*this);
                        return old;
                    }
                }; 
                ViewToVectorOfPointers(std::vector<Entity*>& vec): vec(vec) {}
                auto begin() const { return PointersViewIterator { vec.begin() }; }
                auto end() const { return PointersViewIterator { vec.end() }; }
                std::size_t size() const { return vec.size(); }
            };
            // end of inner class vector_of_pointers_view
            //-------------------------------------------------
            return ViewToVectorOfPointers{
                 m_entitiesByType[typeid(ActualT)] };
        }
   
        // Additional auxiliary functions:

        // (Cell::C1) 
        std::size_t numEntities() const noexcept {
            return m_allEntities.size();
        }

        // (Cell::C2) 
        // NOTE: complexity required to be O(1)
        template<typename ActualT> requires concrete_derived_or_same<ActualT, Entity>
        std::size_t numEntities() const noexcept { 
            // do not create an entry just because of a call to size
            auto itr = m_entitiesByType.find(typeid(ActualT));
            if(itr == m_entitiesByType.end()) return 0;
            return itr->second.size();
        }

        // Iterators begin and end:
        // The Cell would have begin and end iterators for iterating over all pointers to Entity. Retrieved pointers are non-const, i.e. the user can call non-const methods on the retrieved entities. There is no defined order for the iteration. Iteration itself shall not create any copies.
        auto begin() const noexcept { 
            return m_allEntities.begin(); 
        }

        auto end() const noexcept { 
            return m_allEntities.end(); 
        }
    };

    // initializes the cells
    void initializeCells(){
        Latitude bottomLat(-90);
        Latitude topLat(0);
        for(std::size_t curRow = 0; curRow < num_rows; ++curRow)
        {
            topLat = bottomLat + m_cellHeightDeg;
            // Calculate average perimeter
            auto topPerim = CoordinatesMath::perimeterOnLatitude(topLat);
            auto botPerim = CoordinatesMath::perimeterOnLatitude(bottomLat);
            auto avgPerim = (topPerim + botPerim) / 2;
            // get num of cells in row
            int numRowCells = std::ceil(avgPerim / m_cellHeightMeters);
            m_numCells += numRowCells;
            // initialize row struct
            GridRow newRow(Longitude(360 /((double) numRowCells)), std::vector<Cell*>());
            // We add all the cells for our iterator to use
            for(int i = 0; i <  numRowCells; ++i){
                // Add the actual cell to our linked-list
                Cell &p_cell = m_cells.emplace_back(Cell());
                // Add the poiner for the row
                newRow.rowCells.push_back(&p_cell);
            }
            // move the row struct to our vector
            m_gridCells.emplace_back(std::move(newRow));

            // go to next bottom lat
            bottomLat += m_cellHeightDeg;
        }
    }


    // returns the column for row row at longitude lon
    // cell holds: [Left, Right), [Top, Bottom)
    std::size_t rowAndLonToCol(std::size_t row, Longitude lon) const{
        // if we are at the rightmost edge we need the last cell
        if(lon == Longitude(180)) return 0;
        auto colCellWidth = m_gridCells.at(row).cellWidthDeg;
        return std::floor((CoordinatesMath::wrap180(static_cast<double>(lon)) + 180) / static_cast<double>(colCellWidth));
    }

    // Converts latitude to row
    // cell holds: [Left, Right), [Top, Bottom)
    std::size_t latToRow(Latitude lat) const{
        //if we are at the north pole we need to return the last row
        if(lat == Latitude(90)) return m_gridCells.size() - 1;
        // need to align to 180 and divide by cell height to get the row
        return std::floor((CoordinatesMath::wrap90(static_cast<double>(lat)) + 90) / static_cast<double>(m_cellHeightDeg));
    }

    // Non const return of cell at coordinates
    Cell* getCellAtInner(Coordinates c) const{
        // We just get the row and col and return the value
        auto row = latToRow(c.latitude());
        auto col = rowAndLonToCol(row, c.longitude());        
        return m_gridCells.at(row).rowCells.at(col);
    }

    // Truncates a coordinate to cell
    double truncateCoordinates(double coordinate, double precision) const
    {           
            return std::floor(coordinate / precision) * precision;
    }

    // truncates coordinates to representing coordinates of containing cell
    Coordinates truncateCoordinatesToCell(const Coordinates &coords) const
    {
        auto row = latToRow(coords.latitude());
        auto col = rowAndLonToCol(row, coords.longitude());
        auto rowWidth = m_gridCells.at(row).cellWidthDeg;
        Latitude lat =  Latitude(((double) row + 0.5)*  static_cast<double>(m_cellHeightDeg) - 90) ;
        Longitude lon = Longitude(((double) col + 0.5) * static_cast<double>(rowWidth) - 180);
        return Coordinates(lon, lat);

    }

    // Checks if a cell is contained within a given radius from a point
    bool isCellInRadius(const Coordinates &coords, const Coordinates &center, Meters radius) const{
        auto rowWidth = m_gridCells.at(latToRow(coords.latitude())).cellWidthDeg;
        Longitude halfWidth(static_cast<double>(rowWidth) / 2);
        Latitude halfHeight(static_cast<double>(m_cellHeightDeg) / 2); 
        // We just check the if the distance from the center to any of the edges is smaller or equal to radius
        Coordinates bottomLeft = Coordinates(coords.longitude() - halfWidth, coords.latitude() - halfHeight);
        Coordinates bottomRight = Coordinates(coords.longitude() + halfWidth, coords.latitude() - halfHeight);
        Coordinates topLeft = Coordinates(coords.longitude() - halfWidth, coords.latitude() + halfHeight);
        Coordinates topRight = Coordinates(coords.longitude() + halfWidth, coords.latitude() + halfHeight);
        // check if any edge is close enough
        if(CoordinatesMath::distanceFromSegment(center, bottomLeft, topLeft) <= radius) return true;
        if(CoordinatesMath::distanceFromSegment(center, bottomLeft, bottomRight) <= radius) return true;
        if(CoordinatesMath::distanceFromSegment(center, bottomRight, topRight) <= radius) return true;
        if(CoordinatesMath::distanceFromSegment(center, topRight, topLeft) <= radius) return true;
        return false;
    }

    // Aligns coordinates to longitude -180,180 and latitude -90 , 90
    Coordinates alignCoordinates(const Coordinates &coords) const{

        double newLon = CoordinatesMath::wrap180(static_cast<double>(coords.longitude()));
        double newLat = CoordinatesMath::wrap90(static_cast<double>(coords.latitude()));
        return Coordinates(Longitude(newLon), Latitude(newLat));
    }
    
    // adds all the cells in the row corresponding to nextCoordsRow which are at most radius away from center
    void getCellsInRow(const Coordinates &nextCoords, const Coordinates &center,
        Meters radius, std::unordered_set<const Cell*> &visitedCells ,  std::unordered_set<Coordinates> &visitedCoords , 
        std::vector<const Cell*> &results) const{
        auto nextCoordsRow = nextCoords;
        Longitude rowWidth = m_gridCells.at(latToRow(nextCoordsRow.latitude())).cellWidthDeg;
        // Going forward in row
        while(isCellInRadius(alignCoordinates(nextCoordsRow), center, radius) && !visitedCoords.contains(alignCoordinates(nextCoordsRow))){
            // add the cell because it's in range
            auto curCell = getCellAt(nextCoordsRow);
            if(!visitedCells.contains(curCell)) {
                results.emplace_back(curCell); 
                visitedCells.insert(curCell);
            }
            visitedCoords.insert(alignCoordinates(nextCoordsRow));
            // go to the next cell in row
            nextCoordsRow = Coordinates(Longitude(nextCoordsRow.longitude() + rowWidth), nextCoords.latitude());
        }
        // go backwards in row
        nextCoordsRow = Coordinates(Longitude(nextCoords.longitude() - rowWidth), nextCoords.latitude());
        while(isCellInRadius(alignCoordinates(nextCoordsRow), center, radius) && !visitedCoords.contains(alignCoordinates(nextCoordsRow))){
            // add the cell because it's in range
            auto curCell = getCellAt(nextCoordsRow);
            if(!visitedCells.contains(curCell)) {
                results.emplace_back(curCell);
                visitedCells.insert(curCell);
            }
            visitedCoords.insert(alignCoordinates(nextCoordsRow));
            // go to the previous cell in row
            nextCoordsRow = Coordinates(Longitude(nextCoordsRow.longitude() - rowWidth), nextCoords.latitude()); 
        }
    }

    public:
    // Grid’s Constructors and Assignment:
    Grid() : m_cellHeightDeg((double)180 / num_rows), m_cellHeightMeters(2 * CoordinatesMath::half_earth_hemisphere / num_rows), m_numCells(0)
    {
        initializeCells();
    }

    // Grid’s Destructor:
    // You should decide whether there is a need to implement a user defined destructor.

    // Adding Entities to the Grid:

    // (Grid::A1) 
    template<typename ActualT> requires derived_or_same<ActualT, Entity>
    const Cell* add(Coordinates c, ActualT& e) {
        // just get the cell and to it
        auto cell = getCellAtInner(c);
        cell->add(e);
        return cell;
    }

    // Getting Cells from the Grid:
    // (Grid::B1) 
    const Cell* getCellAt(Coordinates c) const {
        return getCellAtInner(c);
    }
    
    // (Grid::B2) (BONUS: you may choose not to implemnt this function)
    std::vector<const Cell*> getCellsAt(Coordinates center, Meters radius) const {
        // This will hold our results
        std::vector<const Cell*> results;
        std::unordered_set<const Cell*> visitedCells;
        std::unordered_set<Coordinates> visitedCoords;
        auto nextCoords = truncateCoordinatesToCell(center);
        // In this loop we are going upwards from center cell and checking adding all the cells in the correct range
        while(isCellInRadius(alignCoordinates(nextCoords), center, radius) && !visitedCoords.contains(alignCoordinates(nextCoords))){
            // get all the cells in the row
            getCellsInRow(nextCoords, center, radius, visitedCells, visitedCoords,  results);
            // go to row above
            nextCoords = Coordinates(nextCoords.longitude(), nextCoords.latitude() + m_cellHeightDeg);
        }
        // going downards now
        nextCoords = truncateCoordinatesToCell(center);
        nextCoords = Coordinates(nextCoords.longitude(), nextCoords.latitude() - m_cellHeightDeg);
        while(isCellInRadius(alignCoordinates(nextCoords), center, radius) && !visitedCoords.contains(alignCoordinates(nextCoords))){
            // get all the cells in the row
            getCellsInRow(nextCoords, center, radius, visitedCells, visitedCoords , results);
            // go to row below
            nextCoords = Coordinates(nextCoords.longitude(), nextCoords.latitude() - m_cellHeightDeg);
            // Handle direction change
        }
        // If the radius is so smaller that we didn't reach any of the edges of the cell containing the coordinate, only this cell is
        // in the range
        if(results.size() == 0) results.emplace_back(getCellAt(center));
        return results;
    }

    // additional auxiliary functions:

    // (Grid::C1) 
    std::size_t numRows() const noexcept {
        return num_rows;
    } 
    
    // (Grid::C2) 
    std::size_t numCols(Coordinates c) const noexcept {
        return m_gridCells.at(latToRow(c.latitude())).rowCells.size();
    }

    // (Grid::C3) 
    std::size_t numCells() const noexcept {
        return m_numCells; 
    }

    // Set this to true if you are implementing a sparse grid.
    static constexpr bool is_sparse = false;  
 
    // Following would iterate over:
    // 1. Only non empty Cells, if is_sparse==true 
    // 2. All Cells, if is_sparse==false
    auto begin() const noexcept { 

        return m_cells.begin(); 
    }

    auto end() const noexcept { 
        return m_cells.end(); 
    }
};