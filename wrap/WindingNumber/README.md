# Fast Winding Numbers for Soups

https://github.com/alecjacobson/WindingNumber

Implementation of the _ACM SIGGRAPH_ 2018 paper, 

"Fast Winding Numbers for Soups and Clouds" 

Gavin Barill¹, Neil Dickson², Ryan Schmidt³, David I.W. Levin¹, Alec Jacobson¹

¹University of Toronto, ²SideFX, ³Gradient Space


_Note: this implementation is for triangle soups only, not point clouds._

This version does _not_ depend on Intel TBB. Instead it depends on
[libigl](https://github.com/libigl/libigl)'s simpler `igl::parallel_for` (which
uses `std::thread`)

<del>This code, as written, depends on Intel's Threading Building Blocks (TBB) library for parallelism, but it should be fairly easy to change it to use any other means of threading, since it only uses parallel for loops with simple partitioning.</del>

The main class of interest is UT_SolidAngle and its init and computeSolidAngle functions, which you can use by including UT_SolidAngle.h, and whose implementation is mostly in UT_SolidAngle.cpp, using a 4-way bounding volume hierarchy (BVH) implemented in the UT_BVH.h and UT_BVHImpl.h headers.  The rest of the files are mostly various supporting code.  UT_SubtendedAngle, for computing angles subtended by 2D curves, can also be found in UT_SolidAngle.h and UT_SolidAngle.cpp .

An example of very similar code and how to use it to create a geometry operator (SOP) in Houdini can be found in the HDK examples (toolkit/samples/SOP/SOP_WindingNumber) for Houdini 16.5.121 and later.  Query points go in the first input and the mesh geometry goes in the second input.


Create a single header using:

    echo "// This header created by issuing: \`$BASH_COMMAND\` $(echo "" | cat - LICENSE README.md | sed -e "s#^..*#\/\/ &#") $(echo "" | cat - SYS_Types.h SYS_Math.h VM_SSEFunc.h VM_SIMD.h UT_Array.h UT_ArrayImpl.h UT_SmallArray.h UT_FixedVector.h UT_ParallelUtil.h UT_BVH.h UT_BVHImpl.h UT_SolidAngle.h UT_Array.cpp UT_SolidAngle.cpp | sed -e "s/^#.*include  *\".*$//g")"
