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
        std::vector<Entity*> getEntities(PredicateFunc&&) const {
            std::vector<Entity*> filteredEntities;
            std::copy_if(m_allEntities.begin(), m_allEntities.end(), std::back_inserter(filteredEntities), 
            [](Entity *entity){ return PredicateFunc(*entity);} );
            return filteredEntities;
        }

        // (Cell::B1) 
        // Getting specific type of Entities from a Cell:
        // returns all entities of type ActualT which return true for the PredicateFunc.
        template<typename ActualT, typename PredicateFunc> 
        requires concrete_derived_or_same<ActualT, Entity>
        std::vector<ActualT*> getEntities(PredicateFunc&&) const {
            std::vector<Entity*> filteredEntities;
            auto typeItr = m_entitiesByType.find(typeid(ActualT));
            std::copy_if(typeItr.begin(), typeItr.end(), std::back_inserter(filteredEntities), 
            [](Entity *entity){ return PredicateFunc(*entity);} );
            return filteredEntities;
        }

        // (Cell::B2) 
        // same as (Cell::B1) above but with a limit on the number of returned entities (up to limit entities).
        template<typename ActualT, typename PredicateFunc> requires concrete_derived_or_same<ActualT, Entity>
        std::vector<ActualT*> getEntities(PredicateFunc&&, std::size_t limit) const {
            std::vector<Entity*> filteredEntities;
            std::size_t curCopied = 0;
            auto typeItr = m_entitiesByType.find(typeid(ActualT));
            // copy untill full or untill no more to copy
            for(auto entity : typeItr)
            {
                if(PredicateFunc(*entity)){
                    filteredEntities.push_back(entity);
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
    std::size_t rowAndLonToCol(std::size_t row, Longitude lon) const{
        auto colCellWidth = m_gridCells.at(row).cellWidthDeg;
        return std::floor((CoordinatesMath::wrap180(static_cast<double>(lon)) + 180) / static_cast<double>(colCellWidth));
    }

    // Converts latitude to row
    std::size_t latToRow(Latitude lat) const{
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
    

    // Truncates a Longitude coordinate to right
    double truncateLongitudeCoordinate(double coordinate, double precision) 
    {           
            return std::floor(coordinate / precision) * precision;
    }

    // Truncates a Latitude to bottom
    double truncateLatitudeCoordinate(double coordinate, double precision) 
    {           
            return std::floor(coordinate / precision) * precision;
    }

    // truncates coordinates to representing coordinates of containing cell
    Coordinates truncateCoordinatesToCell(const Coordinates &coords)
    {
        double rowWidthOrig = m_gridCells.at(latToRow(coords.latitude())).cellWidthDeg;
        double longitude = truncateLongitudeCoordinate(static_cast<double>(coords.longitude()), rowWidthOrig);
        double latitude = truncateLatitudeCoordinate(static_cast<double>(coords.latitude()),static_cast<double>(m_cellHeightDeg));
        // 180 [deg] and - 180 [deg] are the same cell
        if(std::fabs(longitude + 180) <= Epsilon)
        {
            longitude = 180;
        } 
        // If we are at one of the poles are cells
        if((std::fabs(latitude - 90) <= Epsilon) || (std::fabs(latitude + 90) <= Epsilon))
        {
            longitude = 0;
        }
        return {Longitude{(longitude)},
            Latitude{latitude}};
    }

    // Adds all surrounding cell coords to the stack
    void pushCellNeighborsToStack(std::stack<Coordinates> &stack, const Coordinates &cellCoords, std::unordered_set<Coordinates> visitedCoords){
        
        Longitude longitude = cellCoords.longitude();
        Latitude latitude = cellCoords.latitude();
        double latDbl = static_cast<double>(latitude);
        double rowWidthOrig = m_gridCells.at(latToRow(latitude)).cellWidthDeg;

        // Our poles
        static const Coordinates NorthPoleCoords(Longitude(0), Latitude(90));
        static const Coordinates SouthPoleCoords(Longitude(0), Latitude(-90));

        // If we are at a pole, our neighbors are every cell surrounding it
        if(std::fabs(latDbl - 90) <= Epsilon)
        {
            Latitude neighborLat(90 - static_cast<double>(m_cellHeightDeg));
            for(double lon = -180 ; lon <= 180 ; lon += rowWidthOrig)
            {
            Coordinates neighborCoords(Longitude(lon), neighborLat);
                if(!(visitedCoords.contains(neighborCoords)))
                        stack.push(neighborCoords);
            }
        }

        if(std::fabs(latDbl  + 90) <= Epsilon)
        {
            Latitude neighborLat(-90 + static_cast<double>(m_cellHeightDeg));
            for(double lon = -180 ; lon <= 180 ; lon += rowWidthOrig)
            {
            Coordinates neighborCoords(Longitude(lon), neighborLat);
                if(!(visitedCoords.contains(neighborCoords)))
                        stack.push(neighborCoords);
            }
        }
        
        // If we are not at the poles
        if(latDbl < 90 && latDbl > -90){

            Latitude northLat = Latitude(latitude + m_cellHeightDeg);
            Latitude southLat = Latitude(latitude - m_cellHeightDeg);
            // Note: Longitudes are phase-aligned in the constructor
            Longitude eastLon = Longitude(longitude + Longitude(rowWidthOrig));
            Longitude westLon = Longitude(longitude - Longitude(rowWidthOrig));

            bool northPoleTest = (northLat != Latitude(90));
            bool southPoleTest = (southLat != Latitude(-90));
            // If we are at the poles we need to set the longitude to 0.
            Longitude eastLonNorth = northPoleTest ? eastLon : Longitude(0);
            Longitude westLongNorth = northPoleTest ? westLon : Longitude(0);
            Longitude eastLonSouth = southPoleTest ? eastLon : Longitude(0);
            Longitude westLongSouth = southPoleTest ? westLon : Longitude(0);

            Coordinates northCoords = northPoleTest ? Coordinates(longitude, northLat) : NorthPoleCoords;
            Coordinates southCoords = southPoleTest ? Coordinates(longitude, southLat) : SouthPoleCoords;
            Coordinates eastCoords(eastLon,latitude);
            Coordinates westCoords(westLon, latitude);
            Coordinates northWestCoords(westLongNorth, northLat);
            Coordinates northEastCoords(eastLonNorth, northLat);
            Coordinates southWestCoords(westLongSouth, southLat);
            Coordinates southEastCoords(eastLonSouth, southLat);

            // Push
            if(!(visitedCoords.contains(northCoords)))
                    stack.push(northCoords);
            if(!(visitedCoords.contains(southCoords)))
                    stack.push(southCoords);
            if(!(visitedCoords.contains(eastCoords)))
                    stack.push(eastCoords);
            if(!(visitedCoords.contains(westCoords)))
                    stack.push(westCoords);
            if(!(visitedCoords.contains(northWestCoords)))
                    stack.push(northWestCoords);
            if(!(visitedCoords.contains(northEastCoords)))
                    stack.push(northEastCoords);
            if(!(visitedCoords.contains(southWestCoords)))
                    stack.push(southWestCoords);
            if(!(visitedCoords.contains(southEastCoords)))
                    stack.push(southEastCoords);
        }
    }

    // Checks if a cell is contained within a given radius from a point
    bool isCellInRadius(const Coordinates &coords, const Coordinates &center, Meters radius){
        auto rowWidth = m_gridCells.at(latToRow(coords.latitude())).cellWidthDeg; 
        // We just check the if the distance from the center to any of the edges is smaller or equal to radius
        Coordinates bottomRight = Coordinates(coords.longitude() + rowWidth, coords.latitude());
        Coordinates topLeft = Coordinates(coords.longitude(), coords.latitude() + m_cellHeightDeg);
        Coordinates topRight = Coordinates(coords.longitude() + rowWidth, coords.latitude() + m_cellHeightDeg);

        // check if any edge is close enough
        if(CoordinatesMath::distanceFromSegment(center, coords, topLeft) <= radius) return true;
        if(CoordinatesMath::distanceFromSegment(center, coords, bottomRight) <= radius) return true;
        if(CoordinatesMath::distanceFromSegment(center, bottomRight, topRight) <= radius) return true;
        if(CoordinatesMath::distanceFromSegment(center, topRight, topLeft) <= radius) return true;
        return false;
    }

    public:
    // Grid’s Constructors and Assignment:
    Grid() : m_cellHeightDeg((double)180 / num_rows), m_cellHeightMeters((4 * CoordinatesMath::half_earth_hemisphere) / 180), m_numCells(0)
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
        // This is the stack for traversing the cells via coordinates
        std::stack<Coordinates> traversalStack;
        // We dont want to re-visit coordinates
        std::unordered_set<Coordinates> visitedCoords;
        // insert the intial coordinates
        Coordinates initialCoordinates = truncateCoordinatesToCell(center);
        traversalStack.push(initialCoordinates);
        // We traverse while the stack isnt empty
        do{
            // get first coordinates
            auto curCoords = traversalStack.top();
            traversalStack.pop();
            // Ignore entities we already checked
            if(visitedCoords.contains(curCoords)) continue;
            visitedCoords.insert(curCoords);
            // If the cell isnt in the radius we stop, otherwise we continue to the next cells
            if(!isCellInRadius(curCoords, center, radius)) continue;
            // otherwise we insert the cell and check the neighbors
            results.emplace_back(getCellAt(curCoords));
            pushCellNeighborsToStack(visitedCoords, curCoords, visitedCoords);

        }while(!traversalStack.empty());
        // We then traverse in a BFS like fashion on the cells, we stop when the stack is empty
        return std::vector<const Cell*>{};
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
