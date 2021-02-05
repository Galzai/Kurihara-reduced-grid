#pragma once

#include "CoordinatesMath.h"
#include "GISNamedTypes.h"
#include <typeinfo>
#include <typeindex>
#include <concepts>
#include <ranges>
#include <unordered_map>

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
    std::vector<Cell*> m_cellRefs;

    // Simple struct for compact representation of row
    struct GridRow{
        Longitude cellWidthDeg;
        std::vector<Cell> rowCells;
        GridRow(Longitude cellWidthVal, std::vector<Cell> rowCellsVal) : cellWidthDeg(cellWidthVal), rowCells(rowCellsVal){};
    };
    
    // private inner class Cell
    class Cell {
        mutable std::unordered_map<std::type_index, std::vector<Entity*>> m_entitiesByType;
        mutable std::vector<Entity*> m_allEntities;

    public:
        // Add an entity of concrete type to cell by type and to all
        template<typename ActualT> requires derived_or_same<ActualT, Entity>
        void add(ActualT& e) {
            auto itr = m_entitiesByType.find(typeid(ActualT));
            // If theres no type we need to create an entry
            if(itr == m_entitiesByType.end())
            {
                m_entitiesByType.emplace(typeid(ActualT) ,std::vector<Entity*>()); 
            }
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
            // 'm_entitiesByType' is mutable, empty vector might be created below
            auto itr = m_entitiesByType.find(typeid(ActualT));
            // No types means return an empty vector
            if(itr == m_entitiesByType.end())
            {
                m_entitiesByType.emplace(typeid(ActualT) ,std::vector<Entity*>()); 
            }
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
            GridRow newRow(Longitude(360 / numRowCells), std::vector<Cell>(numRowCells, Cell()));
            // We add all the cells for our iterator to use
            for(auto &cell: newRow.rowCells){
                m_cellRefs.emplace_back(&cell);
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
        return &(*(m_gridCells.at(row).rowCells.begin() + col));
    }
    
public:
    // Grid’s Constructors and Assignment:
    Grid() : m_cellHeightDeg(180 / num_rows), m_cellHeightMeters((4 * CoordinatesMath::half_earth_hemisphere) / 180), m_numCells(0)
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
        // TODO: Implement
        (void) center; 
        (void) radius; 
        return std::vector<const Cell*>{};
    }

    // additional auxiliary functions:

    // (Grid::C1) 
    std::size_t numRows() const noexcept {
        return num_rows;
    } 
    
    // (Grid::C2) 
    std::size_t numCols(Coordinates c) const noexcept {
        return m_gridCells.at(latToRow(c.latitude()));
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

        return m_cellRefs.begin(); 
    }

    auto end() const noexcept { 
        return m_cellRefs.end(); 
    }
};
