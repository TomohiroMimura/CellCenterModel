# CellCenterModel

This is a program for a cell-centered model that represents the three-dimensional morphogenesis of a monolayer.

## videos

### In-plane deformation

| with cell rearrangement | without cell rearrangement |
|----------|------|
|[![](https://img.youtube.com/vi/QdJE2k9OQ_Q/maxresdefault.jpg)](https://www.youtube.com/watch?v=QdJE2k9OQ_Q) | [![](https://img.youtube.com/vi/6llLcvThxZA/maxresdefault.jpg)](https://www.youtube.com/watch?v=6llLcvThxZA) |

### out-of-plane deformation

| uniaxial cell division | random cell division |
|----------|------|
| [![](https://img.youtube.com/vi/_A38DCdr_u4/maxresdefault.jpg)](https://www.youtube.com/watch?v=_A38DCdr_u4)|[![](https://img.youtube.com/vi/3tw8Ug5oa7E/maxresdefault.jpg)](https://www.youtube.com/watch?v=3tw8Ug5oa7E) |

### apical constriction

| circular pattern | circumferential pattern|
|--------|---------|
| [![](https://img.youtube.com/vi/AcwwKy-2Oyw/maxresdefault.jpg)](https://www.youtube.com/watch?v=AcwwKy-2Oyw)| [![](https://img.youtube.com/vi/W1OMQtU5rZo/maxresdefault.jpg)](https://www.youtube.com/watch?v=W1OMQtU5rZo)|

## how to make and run ##
Parameter files (parameters.hpp) for the three patterns are included in each directory.

When executing, place one of them in the same hierarchy as the makefile.
Then execute the following code.

    make
    ./a.out




Any published work which utilizes this code shall include the following reference:

Tomohiro Mimura and Yasuhiro Inoue "Cell-center-based model for simulating three-dimensional monolayer tissue deformation", Journal of Theoretical Biology (2023) 
```
@article{MIMURA2023111560,
    title = {Cell-center-based model for simulating three-dimensional monolayer tissue deformation},
    journal = {Journal of Theoretical Biology},
    volume = {571},
    pages = {111560},
    year = {2023},
    issn = {0022-5193},
    author = {Tomohiro Mimura and Yasuhiro Inoue}
}
```
