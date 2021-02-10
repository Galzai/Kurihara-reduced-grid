#include "gtest/gtest.h"
#include "Grid.h"
#include "CoordinatesMath.h"
#include <cassert>


/****************************** Cell Tests ******************************/
// This tests getEntitiesView
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
    Grid<A, 1200> grid;
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

// Checks if getEntities works with abstract type
TEST(GISEx4, getEntitiesPredicateForAbstract){
        struct A {
        virtual ~A() {}
        virtual bool foo() const = 0;
    };
    struct B : A { bool foo() const override {return true;} };
    struct C : A { bool foo() const override {return true;} };
    struct D : A { bool foo() const override {return false;} };
    Grid<A, 400> grid;
    B b1, b2;
    C c;
    D d;

    Coordinates coord{Longitude{45}, Latitude{89.009}};
    const auto p_cell = grid.add(coord, b1);
    grid.add(coord, b2);
    grid.add(coord, c);
    grid.add(coord, d);
    auto predicate = [](const A& en){
        return en.foo();
    };
    EXPECT_EQ(p_cell->getEntities(predicate).size(), (std::size_t) 3);
} 

// Checks if getEntities works with concrete type with predicate signature bool pred(const Entity&);
TEST(GISEx4, getEntitiesAbstractPredicateForConcrete){
        struct A {
        virtual ~A() {}
        virtual bool foo() const = 0;
    };
    struct B : A {B(int val) :x(val){}; int x; bool foo() const override {return x == 1;} };
    struct C : A { bool foo() const override {return true;} };
    struct D : A { bool foo() const override {return false;} };
    Grid<A, 5> grid;
    B b1(1), b2(2) ,b3(1);
    C c;
    D d;

    Coordinates coord{Longitude{45}, Latitude{89.009}};
    const auto p_cell = grid.add(coord, b1);
    grid.add(coord, b2);
    grid.add(coord, b3);
    grid.add(coord, c);
    grid.add(coord, d);
    auto predicate = [](const A& en){
        return en.foo();
    };
    EXPECT_EQ(p_cell->getEntities<B>(predicate).size(), (std::size_t) 2);
} 

// Checks if getEntities works with concrete type with predicate signature bool pred(const ActualT&);
TEST(GISEx4, getEntitiesConcretePredicateForConcrete){
        struct A {
        virtual ~A() {}
        virtual bool foo() const = 0;
    };
    struct B : A {B(int val) :x(val){}; int x; bool foo() const override {return x == 1;} };
    struct C : A { bool foo() const override {return true;} };
    struct D : A { bool foo() const override {return false;} };
    Grid<A, 150> grid;
    B b1(1), b2(2) ,b3(1);
    C c;
    D d;

    Coordinates coord{Longitude{45}, Latitude{89.009}};
    const auto p_cell = grid.add(coord, b1);
    grid.add(coord, b2);
    grid.add(coord, b3);
    grid.add(coord, c);
    grid.add(coord, d);
    auto predicate = [](const B& en){
        return en.foo();
    };
    EXPECT_EQ(p_cell->getEntities<B>(predicate).size(), (std::size_t) 2);
}


// Checks if getEntities works with concrete type with predicate signature bool pred(const ActualT&) and limit;
TEST(GISEx4, getEntitiesConcretePredicateForConcreteAndLimit){
        struct A {
        virtual ~A() {}
        virtual bool foo() const = 0;
    };
    struct B : A {B(int val) :x(val){}; int x; bool foo() const override {return x == 1;} };
    struct C : A { bool foo() const override {return true;} };
    struct D : A { bool foo() const override {return false;} };
    Grid<A, 130> grid;
    B b1(1), b2(1) ,b3(1);
    C c;
    D d;

    Coordinates coord{Longitude{45}, Latitude{89.009}};
    const auto p_cell = grid.add(coord, b1);
    grid.add(coord, b2);
    grid.add(coord, b3);
    grid.add(coord, c);
    grid.add(coord, d);
    auto predicate = [](const B& en){
        return en.foo();
    };
    EXPECT_EQ(p_cell->getEntities<B>(predicate, 1).size(), (std::size_t) 1);
}

// Checks if getEntities works with concrete type with predicate signature bool pred(const Entity&) and limit;
TEST(GISEx4, getEntitiesAbstractPredicateForConcreteAndLimit){
        struct A {
        virtual ~A() {}
        virtual bool foo() const = 0;
    };
    struct B : A {B(int val) :x(val){}; int x; bool foo() const override {return x == 1;} };
    struct C : A { bool foo() const override {return true;} };
    struct D : A { bool foo() const override {return false;} };
    Grid<A, 158> grid;
    B b1(1), b2(1) ,b3(1);
    C c;
    D d;

    Coordinates coord{Longitude{45}, Latitude{89.009}};
    const auto p_cell = grid.add(coord, b1);
    grid.add(coord, b2);
    grid.add(coord, b3);
    grid.add(coord, c);
    grid.add(coord, d);
    auto predicate = [](const A& en){
        return en.foo();
    };
    EXPECT_EQ(p_cell->getEntities<B>(predicate, 1).size(), (std::size_t) 1);
} 

/****************************** Grid Tests ******************************/
// Checks if numRows works as expected
TEST(GISEx4, checkNumRows){
        struct A {
        virtual ~A() {}
        virtual bool foo() const = 0;
    };
    Grid<A, 10> grid;
    EXPECT_EQ(grid.numRows(), (std::size_t) 10);
} 
// Checks if numRows works as expected
// when setting 4 we expect to have less cols near the poles thus having more cells near the equator
TEST(GISEx4, checkNumCells){
        struct A {
        virtual ~A() {}
        virtual bool foo() const = 0;
    };
    Grid<A, 4> grid;
    EXPECT_EQ(grid.numCells(), (std::size_t) 20);
} 

// when setting 4 rows we expect to have less cols near the poles and more near the equator
// expecting ~3 columns near the poles and in each one near the equator ~7, longitude should have no effect
TEST(GISEx4, checkNumCols){
        struct A {
        virtual ~A() {}
        virtual bool foo() const = 0;
    };
    Grid<A, 4> grid;
    Coordinates southPole(Longitude(12), Latitude(-90));
    Coordinates northPole(Longitude(15), Latitude(90));
    Coordinates nearEquator1(Longitude(180), Latitude(-20));
    Coordinates nearEquator2(Longitude(-180), Latitude(20));
    EXPECT_EQ(grid.numCols(northPole), (std::size_t) 3);
    EXPECT_EQ(grid.numCols(southPole), (std::size_t) 3);
    EXPECT_EQ(grid.numCols(nearEquator1), (std::size_t) 7);
    EXPECT_EQ(grid.numCols(nearEquator2), (std::size_t) 7);
} 

// Checks the cell iterator is of expected type and validates all the entries are correct
TEST(GISEx4, checkGridIterator){
        struct A {
        virtual ~A() {}
        virtual void foo() const = 0;
    };
    struct B : A { void foo() const override {} };
    B b1;
    Grid<A, 1> grid;
    auto &firstCell = (*grid.begin());
    Coordinates southPole(Longitude(12), Latitude(-90));
    grid.add(southPole, b1);
    EXPECT_EQ(firstCell.numEntities(), (std::size_t) 1);
}

/****************************** Bonus Tests ******************************/
// Checks entire map for cells (radius is larger than earth)
TEST(GISEx4, checkGetEntitiesRadiusEntireMap){
        struct A {
        virtual ~A() {}
        virtual bool foo() const = 0;
    };
    Grid<A, 400> grid;
    Coordinates c(Longitude(180), Latitude(-20));
    EXPECT_EQ(grid.getCellsAt(c, CoordinatesMath::earth_radius * 100).size(), grid.numCells());
}

// Checks entities within 0 radius (should only return one cell)
TEST(GISEx4, checkGetEntitiesZeroRadius){
        struct A {
        virtual ~A() {}
        virtual bool foo() const = 0;
    };
    Grid<A, 12> grid;
    Coordinates c(Longitude(175), Latitude(-20));
    EXPECT_EQ(grid.getCellsAt(c, Meters(0)).size(), (std::size_t) 1);
} 

// Checks entities within 0 radius on border (should return 2)
TEST(GISEx4, checkGetEntitiesRadiusBorder){
        struct A {
        virtual ~A() {}
        virtual bool foo() const = 0;
    };
    Grid<A, 12> grid;
    Coordinates c(Longitude(180), Latitude(-20));
    EXPECT_EQ(grid.getCellsAt(c, Meters(0)).size(), (std::size_t) 2);
} 

// Checks entities within same radius on poles - all sizes should be identical
TEST(GISEx4, checkGetEntitiesRadiusPoles){
        struct A {
        virtual ~A() {}
        virtual bool foo() const = 0;
    };
    Grid<A, 20> grid;

    Coordinates c1(Longitude(12), Latitude(90));
    Coordinates c2(Longitude(92), Latitude(-90));
    Coordinates c3(Longitude(40), Latitude(90));
    Coordinates c4(Longitude(50), Latitude(-90));

    int size1 = grid.getCellsAt(c1, Meters(2000000)).size();
    int size2 = grid.getCellsAt(c2, Meters(2000000)).size();
    int size3 = grid.getCellsAt(c3, Meters(2000000)).size();
    int size4 = grid.getCellsAt(c4, Meters(2000000)).size();

    EXPECT_EQ(size1, size2);
    EXPECT_EQ(size2, size3);
    EXPECT_EQ(size3, size4);
} 