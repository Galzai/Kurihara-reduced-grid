#include "gtest/gtest.h"
#include "Grid.h"
#include "CoordinatesMath.h"
#include <cassert>

TEST(GISEx4, DocTest) {
    struct A {
        virtual ~A() {}
        virtual void foo() const = 0;
    };
    struct B : A { void foo() const override {} };
    struct C : A { void foo() const override {} };
    struct D : C { void foodie() const {} };

    Coordinates coord{Longitude{20}, Latitude{30}};
    Grid<A, 1917> grid;
    const auto cell_ptr = grid.getCellAt(coord); 
    auto viewD = cell_ptr->getEntitiesView<D>(); // O(1), always, even for size == N
    EXPECT_EQ(viewD.size(), 0ul); // assume passed successfully

    // adding a D object to the Grid, such that it should be added to the cell at coord  
    D dees[4];
    const auto cell_ptr2 = grid.add(coord, dees[0]);
    EXPECT_EQ(viewD.size(), 1ul); // the new object shall be in the view
    EXPECT_EQ(cell_ptr2, cell_ptr);

    // adding 3 D objects to the Grid, such that it should be added to the cell at coord  
    for(int i=1; i<4; ++i) {
        EXPECT_EQ(grid.add(coord, dees[i]), cell_ptr);
    }
    
    EXPECT_EQ(viewD.size(), 4ul); // the new objects shall be in the view

    int i=0;
    for(D* pd: viewD) { // correct type
        pd->foodie();
        EXPECT_EQ(pd, &dees[i]);
        ++i;
    }
}

// Trying to add an entity to cell and validating it was inserted to correct cell
TEST(GISEx4, addToCell){
    struct A {
        virtual ~A() {}
        virtual void foo() const = 0;
    };
    struct B : A { void foo() const override {} };
    Coordinates coord{Longitude{361}, Latitude{1}};
    Grid<A, 2310> grid;
    B b;
    const auto p_cell2 = grid.getCellAt(coord); 
    EXPECT_EQ(p_cell2->numEntities(), (std::size_t)0);
    const auto p_cell = grid.add(coord, b);
    EXPECT_EQ(p_cell, p_cell2);
    EXPECT_EQ(p_cell->numEntities(), (std::size_t)1);
}

// Checks the cell iterator is of expected type and validates all the entries are correct
TEST(GISEx4, checkCellIterator){
        struct A {
        virtual ~A() {}
        virtual void foo() const = 0;
    };
    struct B : A { void foo() const override {} };
    Coordinates coord{Longitude{90}, Latitude{0}};
    B b;
    Grid<A, 300> grid;
    const auto p_cell = grid.add(coord, b);
    EXPECT_EQ(typeid(*(p_cell->begin())), typeid(A*));
}

// Checks the cell iterator is of expected type and validates all the entries are correct
TEST(GISEx4, checkConcreteNumEntities){
        struct A {
        virtual ~A() {}
        virtual void foo() const = 0;
    };
    struct B : A { void foo() const override {} };
    struct C : A { void foo() const override {} };
    Coordinates coord{Longitude{90}, Latitude{0}};
    B b1, b2;
    C c;
    Grid<A, 300> grid;
    const auto p_cell = grid.add(coord, b1);
    grid.add(coord, b2);
    grid.add(coord, c);
    EXPECT_EQ(p_cell->numEntities<C>(), (std::size_t)1);
    EXPECT_EQ(p_cell->numEntities<B>(), (std::size_t)2);
}