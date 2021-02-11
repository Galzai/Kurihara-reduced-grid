# Grid (Advanced topics in programming Ex4)
This an implementation of a grid represention of the earth simliar to Kurihara reduced grid, for Tel-Aviv universitie's advanced topics in programming course. The grid can contains entities of any types and their derivatives.

## Grid
The grid is implemented as non sparse.

### Constructor
The grid contains only one no-arg constructor, the constructor takes the template size and type arguments, sets the "height" of the rows (in degrees) and calls the `initializeCells` function which:  
* Starts from the south pole -(90 deg which is 0 index row) and ends at the the north pole (90 deg).  
* Calculats the width of each cell in row (by taking the avg perimeter of bottom and top of cell and dividing by the height of the grid's rows in meters).
* Creates all the cells of current row in a linked list of all the cells in the grid
* Adds a pointer to the cell in the cells row and column index in a a vector of row structs (each row struct contains the vector of pointers and the width of the row in degrees).    
* The method stated above allows the cell iterator to be directly used as the linked list`s overhead cleanly and without clutter.

### getCellAt
In order to find the cell, we align each coordinate to 0 - 180 deg for latitude and 0
- 360 deg for longitude, we then find the row by dividing the by cell height, find the row's height from our data structure, and then find the column similarly, while taking into account the north pole and 180 longitude.  

### getCellAt - with radius (Bonus) 
* We first find a representive for coordinate for the cell of the given coordinate (center of the cell). 
* We then starting from the first cell move in both direction in the row, adding all the cells in radius range from center. 
* We keep advancing the row and repeating the previous step until the row is out of range. 
* We then move in the reverse direction from the initial cell repeat the process. 
* We keep track of what coordinates we already encountered to prevent circling back.

### numRows
Simply returns the num_rows template argument

### numCols
Finds the row of the coordinate and returs the size of it

### numCells
returns the size of the cell linked list

### iterators
Returns the linked list of cells cbegin and cend iterators accordingly - this allows for a very simple clean implementation that uses built in iterators without too much overhead (only hold pointers in addition to cells) and no need to access multiple different vectors to output iterator.  

## Cell
The cells implements no constructors besides the default constructor.
for simplicity the each cell holds two members, an unordered map from type index to a vector of pointers to the cell`s entities of the type and a vector of pointers to all the entities in the cell, this allows for a simple and clean iterator (using the vectors existing iterator) and allows for simple function implementation.

### getEntities
Uses copy_if to check all the entities in our vector of all entity pointers to copy only the ones the return true on the predicate.

### getEntities - of type ActualT
Simply iterates over the entities of type actualT and adds them to the output if the predicate returns true.

### getEntities - of type ActualT with limit
same as above only we also end if we reached the limit or found all entities (whatever happens first). 

### getEntitiesView
Implemented using the example from the class (O(1) time complexity).

## numEntities
If by type returns the size of the vector in mapping of type, without type returns the size of the vector of all pointers.

### iterators
Simply returns the cell's all entities vector iterator, this allows for a simple and clean implementation with a using the containers already existing iterator with O(1) additional space (and assumingly more efficient then implementing an iterator which stitches together the map).

## Tests
* DocTest - The cell view test given by the instructors. 
* addToCell -validating entities are correctly inserted to cells.
* checkCellIterator - validates the cell iterator derefrence returns the correct type.
* checkConcreteNumEntities - validates numEntities in cell of concrete type.
* getEntitiesPredicateForAbstract - checks if getEntities with abstract type works
* getEntitiesAbstractPredicateForConcrete - Checks if getEntities works with concrete type with predicate signature bool pred(const Entity&).  
* getEntitiesConcretePredicateForConcrete - Checks if getEntities works with concrete type with predicate signature bool pred(const ActualT&).  
* getEntitiesConcretePredicateForConcreteAndLimit - Checks if getEntities works with concrete type with predicate signature bool pred(const ActualT&) and limit.  
* getEntitiesAbstractPredicateForConcreteAndLimit - Checks if getEntities works with concrete type with predicate signature bool pred(const Entity&) and limit.
* checkNumRows - checks if numRows works as expected.
* checkNumCells - checks num cell works as expected.
* checkNumCols - checks if numCols works as expected, expecting more cols near the equator and less near the poles while being sized symetrically.
* checkGridIterator - checks the grid iterator returns a correct refrence to the cells.  
* checkGetEntitiesRadiusEntireMap -  checks that entering a radius larger than the entire map in getCellsAt returns all tne cells.  
* checkGetEntitiesZeroRadius - checks the setting 0 radius on coordinate that isnt on edge returns only 1 cell in getCellsAt.  
* checkGetEntitiesRadiusBorder - checks that if the radius is 0 and the coordinate is on an edge we get 2 cells in getCellsAt.
* checkGetEntitiesRadiusPoles - checks that at the poles at any longitude, if the radius is identical we get the same number of cells in getCellsAt.  




