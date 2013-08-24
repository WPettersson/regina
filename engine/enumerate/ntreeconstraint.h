
/**
 * A base class for additional linear constraints that we can add to the
 * tableaux of normal surface matching equations.  This is used with
 * NTreeEnumeration, NTreeSingleSoln and related algorithms for enumerating and
 * locating normal surfaces in a 3-manifold triangulation.  See the
 * LPInitialTableaux class notes for details on how these constraints
 * interact with the tableaux of matching equations.
 *
 * The linear constraints may be equalities or inequalities, and there
 * may be more than one such constraint.  If all constraints are
 * homogeneous equalities, the class should derive from LPConstraintSubspace
 * instead (not this base class).
 *
 * This base class provides no functionality.  For documentation's sake
 * only, the notes here describe the functionality that any subclass
 * \e must implement.  We note again that LPConstraintBase does not
 * provide any implementations at all, and subclasses are completely
 * responsible for their own implementations.
 */
class LPConstraintBase {
#ifdef __DOXYGEN
    public:
        enum {
            /**
             * The number of additional linear constraints that we impose.
             * Each constraint will generate one new variable (column)
             * and one new equation (row) in the tableaux.
             */
            nConstraints
        };

        /**
         * Stores the extra coefficients in a single column for the
         * \a nConstraints additional rows that we add to the tableaux
         * to describe the \a nConstraints additional linear equations
         * or inequalities.
         *
         * Subclasses may store these coefficients however they like
         * (in particular, they may optimise for sparse coefficients,
         * binary coefficients, and so on).  They will only ever be
         * accessed through the member functions of this Coefficients class.
         */
        struct Coefficients {
            /**
             * Creates an uninitialised set of coefficients for a single
             * column.  These cofficients must be initialised through a
             * call to addRows() before they can be used.
             */
            Coefficients();

            /**
             * Explicitly fills the final row(s) of the given tableaux matrix 
             * with the coefficients stored in this Coefficients structure.
             * In essence, this routine simply copies this sparse and/or
             * specialised representation of the final row(s) into a
             * more standard dense matrix representation.
             *
             * This routine should only affect the final \a nConstraints
             * entries in the given column of the matrix.  It may assume
             * that these final row(s) have already been initialised to zero.
             *
             * \pre The given matrix has at least \a nConstraints rows
             * and at least \a col + 1 columns.
             * \pre The final \a nConstraints entries in column \a col
             * of the given matrix have already been set to zero.
             *
             * @param m the matrix in which to place these column
             * coefficients.
             * @param col the column of the given matrix in which to
             * place these coefficients.
             */
            void fillFinalRows(LPMatrix& m, unsigned col) const;

            /**
             * Computes the inner product of (i) the final \a nConstraints
             * entries in the given row of the given matrix with (ii) the
             * \a nConstraints column coefficients stored in this data
             * structure.
             *
             * \pre The given matrix has at least \a nConstraints columns
             * and at least \a mRow + 1 rows.
             *
             * @param m the matrix whose row we will use in the inner product.
             * @param mRow the row of the matrix \a m to use in the inner
             * product.
             * @return the resulting portion of the inner product.
             */
            NInteger innerProduct(const LPMatrix& m, unsigned mRow) const;

            /**
             * A variant of innerProduct() that takes into account any
             * adjustments to these linear constraint(s) that are required when
             * this is a quadrilateral column being used to represent an
             * octagon type.
             *
             * The LPData class offers support for octagonal almost normal
             * surfaces, in which exactly one tetrahedron is allowed to have
             * exactly one octagon type.  We represent such an octagon as a
             * \e pair of incompatible quadrilaterals within the same
             * tetrahedron.  See the LPData class notes for details on how
             * this works.
             *
             * In some settings, our extra linear constraints must behave
             * differently in the presence of octagons (i.e., the coefficient
             * of the octagon type is not just the sum of coefficients of the
             * two constituent quadrilateral types).  This routine effectively
             * allows us to adjust the tableaux accordingly.
             *
             * Specifically: this routine computes the inner product of (i) the
             * final \a nConstraints entries in the given row of the given
             * matrix with (ii) the \a nConstraints column coefficients
             * stored in this data structure.  We assume that this column
             * in the underlying tableaux describes one of the two
             * quadrilateral coordinates in some tetrahedron that together
             * form an octagon type, and if necessary we implicitly adjust
             * the coefficients stored in this data structure accordingly.
             *
             * \pre The given matrix has at least \a nConstraints columns
             * and at least \a mRow + 1 rows.
             *
             * \pre This column of the underlying tableaux describes one
             * of the two quadrilateral coordinates that are being
             * combined to form an octagon type within some tetrahedron.
             *
             * @param m the matrix whose row we will use in the inner product.
             * @param mRow the row of the matrix \a m to use in the inner
             * product.
             * @return the resulting portion of the inner product.
             */
            NInteger innerProductOct(const LPMatrix& m, unsigned mRow)
                const;
        };

        /**
         * Explicitly constructs equations for the linear function(s)
         * constrained by this class.  Specifically, this routine takes an
         * array of Coefficients objects (one for each column of the initial
         * tableaux) and fills in the necessary coefficient data.
         *
         * The precise form of the linear function(s) will typically
         * depend upon the underlying triangulation.  For this reason,
         * the triangulation is explicitly passed, along with the
         * permutation that indicates which columns of the initial tableaux
         * correspond to which normal coordinates.
         *
         * More precisely: recall that, for each linear function, the initial
         * tableaux acquires one new variable \a x_i that evaluates this linear
         * function f(x).  This routine must create the corresponding row that
         * sets <tt>f(x) - x_i = 0</tt>.  Thus it must construct the
         * coefficients of f(x) in the columns corresponding to normal
         * coordinates, and it must also set a coefficient of -1 in the
         * column for the corresponding new variable.
         *
         * This function is templated, since in reality we typically
         * pass an array of full tableaux columns (of type
         * LPInitialTableaux::Col), which are larger subclasses of the
         * Coefficients class.  This templating is necessary because the
         * compiler must know how large each column object is in
         * order to correct access each element of the given array.
         *
         * As described in the LPInitialTableaux class notes, it might
         * not be possible to construct the linear functions (since the
         * triangulation might not satisfy the necessary requirements).
         * In this case, this routine should ensure that the linear
         * functions are in fact the zero functions, and should return
         * \c false (but it must still set -1 coefficients for the new
         * variables as described above).  Otherwise (if the linear function
         * were successfully constructed) this routine should return \c true.
         *
         * \pre The template class \a ColClass is a subclass of
         * Coefficients.
         *
         * \pre For all coefficients in the array \a col, the
         * Coefficients substructures have all been initialised with the
         * default constructor and not modified since.
         *
         * @param col the array of columns as stored in the initial
         * tableaux (i.e., the data member LPInitialTableaux::col_).
         * @param columnPerm the corresponding permutation of columns
         * that describes how columns of the tableaux correspond to
         * normal coordinates in the underlying triangulation (i.e., the
         * data member LPInitialTableaux::columnPerm_).
         * @param tri the underlying triangulation.
         * @return \c true if the linear functions were successfully
         * constructed, or \c false if not (in which case they will be
         * replaced with the zero functions instead).
         */
        template <typename ColClass>
        static bool addRows(ColClass* col, const int* columnPerm,
            NTriangulation* tri);

        /**
         * Explicitly constraints each of these linear functions to an
         * equality or inequality in the underlying tableaux.  This will
         * typically consist of a series of calls to LPData::constrainZero()
         * and/or LPData::constrainPositive().
         *
         * The variables for these extra linear functions are stored in
         * columns <tt>numCols - nConstraints</tt>, ..., <tt>numCols - 1</tt>
         * of the given tableaux, and so your calls to LPData::constrainZero()
         * and/or LPData::constrainPositive() should operate on these
         * (and only these) columns.
         * 
         * \pre These column coefficients belong to the initial starting
         * tableaux (LPInitialTableaux) from which the given tableaux is
         * derived.
         *
         * @param lp the tableaux in which to constrain these linear
         * functions.
         * @param numCols the number of columns in the given tableaux.
         */
        static void constrain(LPData<LPConstraintNone>& lp, unsigned numCols);

        /**
         * Ensures that the given normal surface satisfies the extra
         * constraints described by this class.
         *
         * Ideally this test is not based on explicitly recomputing the
         * linear function(s), but instead runs independent tests.
         * For instance, if this class is used to constraint Euler
         * characteristic, then ideally this routine would call
         * s->getEulerCharacteristic() and test the return value of that
         * routine instead.
         *
         * @param s the surface to test.
         * @return \c true if the given surface satisfies these linear
         * constraints, or \c false if it does not.
         */
        static bool verify(const NNormalSurface* s);
#endif
};

/**
 * A subclass of LPConstraintBase used for constraints defined entirely
 * by homogeneous linear equations.
 *
 * Any set of constraints defined entirely by homogeneous linear
 * equations should derive from LPConstraintSubspace, not LPConstraintBase.
 * In other words, any set of constraints derived from LPConstraintSubspace
 * should simply restrict our attention to a vector subspace of the
 * normal surface coordinate system.
 *
 * This class does not provide any additional functionality.  It is
 * merely a convenience to help describe and enforce preconditions.
 */
class LPConstraintSubspace : public LPConstraintBase {
};

/**
 * A do-nothing class that imposes no additional linear constraints on
 * the tableaux of normal surface matching equations.
 */
class LPConstraintNone : public LPConstraintSubspace {
    public:
        enum { nConstraints = 0 };

        struct Coefficients {
            inline Coefficients() {}

            inline void fillFinalRows(LPMatrix& m, unsigned col) const {}

            inline NInteger innerProduct(const LPMatrix&, unsigned) const {
                return NInteger(); // Returns zero.
            }

            inline NInteger innerProductOct(const LPMatrix&, unsigned)
                    const {
                return NInteger(); // Returns zero.
            }
        };

        template <typename ColClass>
        inline static bool addRows(ColClass*, const int*, NTriangulation*) {
            return true;
        }

        inline static void constrain(LPData<LPConstraintNone>&, unsigned) {
        }

        inline static bool verify(const NNormalSurface*) {
            return true;
        }
};

/**
 * A class that constraints the tableaux of normal surface matching equations
 * to ensure that Euler characteristic is strictly positive.
 *
 * There are many ways of writing Euler characteritic as a linear
 * function.  The function constructed here has integer coefficients,
 * but otherwise has no special properties of note.
 *
 * This constraint can work with either normal or almost normal
 * coordinates.  In the case of almost normal coordinates, the function
 * is modified to measure Euler characteristic minus the number of
 * octagons (a technique of Casson, also employed by Jaco and Rubinstein, that
 * is used to ensure we do not have more than two octagons when searching for
 * a normal or almost normal sphere in the 3-sphere recognition algorithm).
 *
 * \pre We are working in standard normal or almost normal coordinates
 * (not quadrilateral or quadrilateral-octagon coordinates).  In
 * particular, the coordinate system passed to the corresponding
 * LPInitialTableaux class constructor must be NNormalSurfaceList::STANDARD.
 */
class LPConstraintEuler : public LPConstraintBase {
    public:
        enum { nConstraints = 1 };

        struct Coefficients {
            int euler;
                /**< The coefficient of the Euler characteristic
                     function for the corresponding column of the matching
                     equation matrix. */

            inline Coefficients() : euler(0) {}

            inline void fillFinalRows(LPMatrix& m, unsigned col) const {
                m.entry(m.rows() - 1, col) = euler;
            }

            inline NInteger innerProduct(const LPMatrix& m,
                    unsigned mRow) const {
                NInteger ans(m.entry(mRow, m.rows() - 1));
                ans *= euler;
                return ans;
            }

            inline NInteger innerProductOct(const LPMatrix& m,
                    unsigned mRow) const {
                // This is called for *two* quad columns (the two quads
                // that combine to give a single octagon).
                //
                // The adjustment in this case is to subtract two from
                // the overall Euler characteristic coefficient for this
                // octagon type (-1 because an octagon has lower Euler
                // characteristic than two quads, and -1 again because
                // we are measuring Euler - #octagons.
                //
                // Happily we can do this by subtracting one from the
                // coefficient in each of the two columns, as
                // implemented below.
                NInteger ans(m.entry(mRow, m.rows() - 1));
                ans *= (euler - 1);
                return ans;
            }
        };

        template <typename ColClass>
        static bool addRows(ColClass* col, const int* columnPerm,
                NTriangulation* tri) {
            int* obj = new int[7 * tri->getNumberOfTetrahedra()];
            unsigned tet, i;
            NPerm4 p;
            for (i = 0; i < 7 * tri->getNumberOfTetrahedra(); ++i)
                obj[i] = 1;
            for (i = 0; i < tri->getNumberOfFaces(); ++i) {
                tet = tri->tetrahedronIndex(
                    tri->getFace(i)->getEmbedding(0).getTetrahedron());
                p = tri->getFace(i)->getEmbedding(0).getVertices();
                --obj[7 * tet + p[0]];
                --obj[7 * tet + p[1]];
                --obj[7 * tet + p[2]];
                --obj[7 * tet + 4];
                --obj[7 * tet + 5];
                --obj[7 * tet + 6];
            }
            for (i = 0; i < tri->getNumberOfEdges(); ++i) {
                tet = tri->tetrahedronIndex(
                    tri->getEdge(i)->getEmbedding(0).getTetrahedron());
                p = tri->getEdge(i)->getEmbedding(0).getVertices();
                ++obj[7 * tet + p[0]];
                ++obj[7 * tet + p[1]];
                ++obj[7 * tet + 4 + vertexSplitMeeting[p[0]][p[1]][0]];
                ++obj[7 * tet + 4 + vertexSplitMeeting[p[0]][p[1]][1]];
            }

            for (i = 0; i < 7 * tri->getNumberOfTetrahedra(); ++i)
                col[i].euler = obj[columnPerm[i]];

            col[7 * tri->getNumberOfTetrahedra()].euler = -1;

            delete[] obj;
            return true;
        }

        inline static void constrain(LPData<LPConstraintEuler>& lp,
                unsigned numCols) {
            lp.constrainPositive(numCols - 1);
        }

        inline static bool verify(const NNormalSurface* s) {
            return (s->getEulerCharacteristic() > 0);
        }
};

/**
 * A class that constraints the tableaux of normal surface matching equations
 * to ensure that normal surfaces in an ideal triangulation are compact
 * (thereby avoiding spun normal surfaces with infinitely many triangles).
 *
 * At present this class can only work with oriented triangulations that have
 * precisely one vertex, which is ideal with torus link.  These
 * constraints are explicitly checked by addRows(), which returns \c false
 * if they are not satisfied.  Moreover, this constraint calls on
 * SnapPea for some calculations: in the unexpected situation where
 * SnapPea retriangulates, the linear function cannot be constructed and
 * addRows() will again return \c false.  You should always test
 * LPInitialTableaux::constraintsBroken() to verify that the linear
 * functions have been constructed correctly.
 *
 * Also, at present this class can only work with quadrilateral normal
 * coordinates (and cannot handle almost normal coordinates at all).
 * This is \e not explicitly checked; instead it appears as a
 * precondition (see below).
 *
 * \pre We are working in quadrilateral normal coordinates.  In particular,
 * the coordinate system passed to the corresponding LPInitialTableaux class
 * must be NNormalSurfaceList::QUAD, and constrainOct() must never be
 * called on any of the corresponding LPData tableaux.
 */
class LPConstraintNonSpun : public LPConstraintSubspace {
    public:
        enum { nConstraints = 2 };

        struct Coefficients {
            int meridian;
            int longitude;

            inline Coefficients() : meridian(0), longitude(0) {}

            inline void fillFinalRows(LPMatrix& m, unsigned col) const {
                m.entry(m.rows() - 2, col) = meridian;
                m.entry(m.rows() - 1, col) = longitude;
            }

            inline NInteger innerProduct(const LPMatrix& m,
                    unsigned mRow) const {
                NInteger ans1(m.entry(mRow, m.rows() - 2));
                ans1 *= meridian;
                NInteger ans2(m.entry(mRow, m.rows() - 1));
                ans2 *= longitude;
                ans1 += ans2;
                return ans1;
            }

            inline NInteger innerProductOct(const LPMatrix& m,
                    unsigned mRow) const {
                // This should never be called, since we never use this
                // constraint with almost normal surfaces.
                // For compilation's sake though, just return the usual
                // inner product.
                return innerProduct(m, mRow);
            }
        };

        template <typename ColClass>
        static bool addRows(ColClass* col, const int* columnPerm,
                NTriangulation* tri) {
            // Regardless of whether the constraints are broken,
            // we need to ensure that the matrix has full rank.
            // Therefore add the coefficients for the two new variables now.
            col[3 * tri->getNumberOfTetrahedra()].meridian = -1;
            col[3 * tri->getNumberOfTetrahedra() + 1].longitude = -1;

            // For the time being we insist on one vertex, which must be
            // ideal with torus link.
            if (tri->getNumberOfVertices() != 1 ||
                    (! tri->getVertex(0)->isIdeal()) ||
                    (! tri->getVertex(0)->isLinkOrientable()) ||
                    tri->getVertex(0)->getLinkEulerCharacteristic() != 0)
                return false;

            // Compute the two slope equations for the torus cusp, if we can.
            NSnapPeaTriangulation snapPea(*tri, false);
            NMatrixInt* coeffs = snapPea.slopeEquations();
            if (! coeffs)
                return false;

            // Check that SnapPea hasn't changed the triangulation on us.
            if (! snapPea.verifyTriangulation(*tri)) {
                delete coeffs;
                return false;
            }

            // All good!  Add the two slope equations as extra rows to
            // our constraint matrix.
            //
            // The coefficients here are differences of terms from
            // SnapPy's get_cusp_equation(), which works in native
            // integers; therefore we will happily convert them back to
            // native integers now.
            for (int i = 0; i < 3 * tri->getNumberOfTetrahedra(); ++i) {
                col[i].meridian = coeffs->entry(0, columnPerm[i]).longValue();
                col[i].longitude = coeffs->entry(1, columnPerm[i]).longValue();
            }

            delete coeffs;
            return true;
        }

        inline static void constrain(LPData<LPConstraintNonSpun>& lp,
                unsigned numCols) {
            lp.constrainZero(numCols - 2);
            lp.constrainZero(numCols - 1);
        }

        inline static bool verify(const NNormalSurface* s) {
            return s->isCompact();
        }
};

/**
 * A base class for additional banning and marking constraints that we
 * can place on tree traversal algorithms.  This is used with
 * NTreeEnumeration, NTreeSingleSoln and related algorithms for
 * enumerating and locating normal surfaces in a 3-manifold triangulation.
 *
 * This class adds constraints of two types:
 *
 * - \e Banning constraints, which ensure that certain normal coordinates
 *   are set to zero;
 *
 * - \e Marking constraints, which are more flexible and can be used in
 *   different ways by different algorithms.
 *
 * All of these constraints operate only on normal coordinates in the
 * underlying tableaux (and in particular not the additional variables
 * introduced by additional linear constraints, as described by
 * LPConstraintBase and its subclasses).
 *
 * Currently marking is used in the following ways:
 *
 * - The NTreeEnumeration algorithm does not use marking at all.
 *
 * - In the NTreeSingleSoln algorithm, marking affects what is considered
 *   a non-trivial surface.  Normally, a non-trivial surface is defined
 *   to be one in which some triangle coordinate is zero.  With marking,
 *   a non-trivial surface is redefined to be one in which some \e unmarked
 *   triangle coordinate is zero.  In other words, marked triangle types
 *   are effectively ignored when determining whether a surface is non-trivial
 *   or not.
 *
 * At present, marking is not used at all for quadrilateral coordinates.
 * Howver, marking is a very new feature, and this concept may be expanded
 * in future versions of Regina.
 *
 * This class does not record disc types in the order of their normal
 * coordinates; instead it records them in the order of their columns in
 * a tableaux for linear programming (as used in LPInitialTableaux).
 * This means that there is a little more work required in setting up
 * the initial lists of banned and marked columns, but then these lists are
 * easy to use on the fly during tree traversal algorithms.
 *
 * This base class provides limited functionality (as documented below).
 * Subclasses \e must implement a constructor (which, like this base
 * class, takes a triangulation and a coordinate system), and must also
 * implement the routine init() which determines which normal coordinates
 * are banned and/or marked.
 */
class BanConstraintBase {
    protected:
        NTriangulation* tri_;
            /**< The triangulation with which we are working. */
        int coords_;
            /**< The normal or almost normal coordinate system in which
                 we are working.  This must be one of NNormalSurfaceList::QUAD,
                 NNormalSurfaceList::STANDARD, NNormalSurfaceList::AN_QUAD_OCT,
                 or NNormalSurfaceList::AN_STANDARD. */
        bool* banned_;
            /**< Indicates which columns of a tableaux correspond to banned
                 disc types.
                 The size of this array is the number of normal coordinates
                 (so we explicitly exclude extra columns that arise from the
                 template parameter LPConstraint. */
        bool* marked_;
            /**< Indicates which columns of a tableaux correspond to marked
                 disc types.
                 The size of this array is the number of normal coordinates
                 (so we explicitly exclude extra columns that arise from the
                 template parameter LPConstraint. */

    protected:
        /**
         * Constructs and initialises the \a banned_ and \a marked_ arrays
         * to be entirely \c false.  The only purpose of passing the
         * triangulation and coordinate system is to determine how many
         * normal coordinates we are dealing with.
         *
         * \warning Before you use this object, the routine init() must be
         * called to fill in the \a banned_ and \a marked_ arrays with the
         * correct data.  Otherwise you will have no banned or marked disc
         * types at all.
         *
         * @param tri the triangulation with which we are working.
         * @param coords the normal or almost normal coordinate system in
         * which we are working.  This must be one of NNormalSurfaceList::QUAD,
         * NNormalSurfaceList::STANDARD, NNormalSurfaceList::AN_QUAD_OCT, or
         * NNormalSurfaceList::AN_STANDARD.
         */
        inline BanConstraintBase(NTriangulation* tri, int coords) :
                tri_(tri), coords_(coords) {
            unsigned nCols = (coords == NNormalSurfaceList::QUAD ||
                coords == NNormalSurfaceList::AN_QUAD_OCT ?
                3 * tri->getNumberOfTetrahedra() :
                7 * tri->getNumberOfTetrahedra());
            banned_ = new bool[nCols];
            marked_ = new bool[nCols];
            std::fill(banned_, banned_ + nCols, false);
            std::fill(marked_, marked_ + nCols, false);
        }

        /**
         * Destroys this object and all associated data.
         */
        inline ~BanConstraintBase() {
            delete[] banned_;
            delete[] marked_;
        }

        /**
         * Enforces all bans described by this class in the given
         * tableaux.  Specifically, for each banned disc type, this
         * routine calls LPData::constrainZero() on the corresponding
         * normal coordinate column.
         *
         * @param lp the tableaux in which to enforce the bans.
         */
        template <typename LPConstraint>
        void enforceBans(LPData<LPConstraint>& lp) const {
            for (unsigned i = 0; i < lp.coordinateColumns(); ++i)
                if (banned_[i])
                    lp.constrainZero(i);
        }

#ifdef __DOXYGEN
        /**
         * Idetifies which disc types to ban and mark, and records the
         * corresponding tableaux columns in the \a banned_ and \a marked_
         * arrays respectively.
         *
         * @param columnPerm the permutation of columns that describes how
         * columns of the tableaux correspond to normal coordinates in
         * the underlying triangulation.  Specifically, this permutation must
         * be the same permutation returned by LPInitialTableaux::columnPerm().
         */
        void init(const int* columnPerm);
#endif
};

/**
 * A do-nothing class that bans no disc types and marks no disc types.
 */
class BanNone : public BanConstraintBase {
    protected:
        /**
         * Constructs and initialises the \a banned_ and \a marked_ arrays
         * to be entirely \c false, as described in the BanConstraintBase
         * superclass constructor.
         *
         * Although one should normally call the routine init() before
         * using this object, for BanNone this is not strictly necessary
         * since there are no disc types to ban or mark.
         *
         * @param tri the triangulation with which we are working.
         * @param coords the normal or almost normal coordinate system in
         * which we are working.  This must be one of NNormalSurfaceList::QUAD,
         * NNormalSurfaceList::STANDARD, NNormalSurfaceList::AN_QUAD_OCT, or
         * NNormalSurfaceList::AN_STANDARD.
         */
        inline BanNone(NTriangulation* tri, int coords) :
                BanConstraintBase(tri, coords) {
        }

        inline void init(const int*) {
        }
};

/**
 * A class that bans normal disc types that meet the boundary of the
 * underlying triangulation.  No disc types are marked at all.
 *
 * \warning This class only works as expected in \e standard normal or
 * almost normal coordinates.  In quadrilateral or quadrilateral-octagon
 * coordinates it will only ban quadrilaterals or octagons that touch
 * the boundary, but it will still allow \e triangles that meet the boundary
 * (since triangle types are not counted in these coordinate systems).
 */
class BanBoundary : public BanConstraintBase {
    protected:
        /**
         * Constructs and initialises the \a banned_ and \a marked_ arrays
         * to be entirely \c false, as described in the BanConstraintBase
         * superclass constructor.
         *
         * \warning Before you use this object, the routine init() must be
         * called to fill in the \a banned_ and \a marked_ arrays with the
         * correct data.  Otherwise you will have no banned or marked disc
         * types at all.
         *
         * @param tri the triangulation with which we are working.
         * @param coords the normal or almost normal coordinate system in
         * which we are working.  This must be one of NNormalSurfaceList::QUAD,
         * NNormalSurfaceList::STANDARD, NNormalSurfaceList::AN_QUAD_OCT, or
         * NNormalSurfaceList::AN_STANDARD.
         */
        inline BanBoundary(NTriangulation* tri, int coords) :
                BanConstraintBase(tri, coords) {
        }

        void init(const int* columnPerm) {
            unsigned n = tri_->getNumberOfTetrahedra();
            unsigned tet, type, i, k;

            bool quadOnly = (coords_ == NNormalSurfaceList::QUAD ||
                coords_ == NNormalSurfaceList::AN_QUAD_OCT);

            // The implementation here is a little inefficient (we repeat tests
            // three or four times over), but this routine is only called at
            // the beginning of the enumeration process so no need to worry.

            // Ban quadrilaterals in tetrahedra that meet the boundary
            // (every such quadrilateral meets a boundary face).
            for (i = 0; i < 3 * n; ++i) {
                if (quadOnly)
                    tet = columnPerm[i] / 3;
                else
                    tet = columnPerm[i] / 7;

                for (k = 0; k < 4; ++k)
                    if (! tri_->getTetrahedron(tet)->adjacentTetrahedron(k)) {
                        banned_[i] = true;
                        break;
                    }
            }

            // Ban triangles in tetrahedra that meet the boundary (but
            // only those triangles that meet the boundary faces).
            if (! quadOnly)
                for (i = 3 * n; i < 7 * n; ++i) {
                    tet = columnPerm[i] / 7;
                    type = columnPerm[i] % 7;

                    for (k = 0; k < 4; ++k)
                        if (k != type &&
                                ! tri_->getTetrahedron(tet)->
                                adjacentTetrahedron(k)) {
                            banned_[i] = true;
                            break;
                        }
                }
        }
};

/**
 * A class that bans and marks disc types associated with torus boundary
 * components.  Here we refer exclusively to real torus boundary
 * components (not ideal vertices with torus cusps).  Specifically:
 *
 * - this class bans any normal triangle or quadrilateral that meets a
 *   torus boundary;
 *
 * - this class marks any normal triangle in the link of a vertex on a
 *   torus boundary.
 *
 * \warning As with BanBoundary, this class only works as expected in
 * \e standard normal or almost normal coordinates.  In quadrilateral or
 * quadrilateral-octagon coordinates it will only ban quadrilaterals or
 * octagons that touch torus boundaries, but it will still allow \e triangles
 * that meet torus boundaries (since triangle types are not counted in these
 * coordinate systems).
 */
class BanTorusBoundary : public BanConstraintBase {
    protected:
        /**
         * Constructs and initialises the \a banned_ and \a marked_ arrays
         * to be entirely \c false, as described in the BanConstraintBase
         * superclass constructor.
         *
         * \warning Before you use this object, the routine init() must be
         * called to fill in the \a banned_ and \a marked_ arrays with the
         * correct data.  Otherwise you will have no banned or marked disc
         * types at all.
         *
         * @param tri the triangulation with which we are working.
         * @param coords the normal or almost normal coordinate system in
         * which we are working.  This must be one of NNormalSurfaceList::QUAD,
         * NNormalSurfaceList::STANDARD, NNormalSurfaceList::AN_QUAD_OCT, or
         * NNormalSurfaceList::AN_STANDARD.
         */
        inline BanTorusBoundary(NTriangulation* tri, int coords) :
                BanConstraintBase(tri, coords) {
        }

        void init(const int* columnPerm) {
            unsigned n = tri_->getNumberOfTetrahedra();
            unsigned tet, type, i, k;

            // Which boundary faces are we banning?
            unsigned nFaces = tri_->getNumberOfFaces();
            bool* banFace = new bool[nFaces];
            std::fill(banFace, banFace + nFaces, false);

            // Which vertex links are we marking triangles around?
            unsigned nVertices = tri_->getNumberOfVertices();
            bool* markVtx = new bool[nVertices];
            std::fill(markVtx, markVtx + nVertices, false);

            NBoundaryComponent* bc;
            for (i = 0; i < tri_->getNumberOfBoundaryComponents(); ++i) {
                bc = tri_->getBoundaryComponent(i);
                if ((! bc->isIdeal()) && bc->isOrientable() &&
                        bc->getEulerCharacteristic() == 0) {
                    // We've found a real torus boundary.
                    for (k = 0; k < bc->getNumberOfFaces(); ++k)
                        banFace[bc->getFace(k)->markedIndex()] = true;
                    for (k = 0; k < bc->getNumberOfVertices(); ++k)
                        markVtx[bc->getVertex(k)->markedIndex()] = true;
                }
            }

            bool quadOnly = (coords_ == NNormalSurfaceList::QUAD ||
                coords_ == NNormalSurfaceList::AN_QUAD_OCT);

            // The implementation here is a little inefficient (we repeat tests
            // three or four times over), but this routine is only called at
            // the beginning of the enumeration process so no need to worry.

            // Ban quadrilaterals that touch torus boundaries.
            for (i = 0; i < 3 * n; ++i) {
                if (quadOnly)
                    tet = columnPerm[i] / 3;
                else
                    tet = columnPerm[i] / 7;

                for (k = 0; k < 4; ++k)
                    if (banFace[tri_->getTetrahedron(tet)->getFace(k)->
                            markedIndex()]) {
                        banned_[i] = true;
                        break;
                    }
            }

            // Ban triangles that touch torus boundaries, and mark all
            // triangles that surround vertices on torus boundaries
            // (even if the triangles do not actually touch the boundary).
            if (! quadOnly)
                for (i = 3 * n; i < 7 * n; ++i) {
                    tet = columnPerm[i] / 7;
                    type = columnPerm[i] % 7;

                    if (markVtx[tri_->getTetrahedron(tet)->getVertex(type)->
                            markedIndex()])
                        marked_[i] = true;

                    for (k = 0; k < 4; ++k)
                        if (k != type &&
                                banFace[tri_->getTetrahedron(tet)->getFace(k)->
                                markedIndex()]) {
                            banned_[i] = true;
                            break;
                        }
                }

            delete[] banFace;
            delete[] markVtx;
        }
};
