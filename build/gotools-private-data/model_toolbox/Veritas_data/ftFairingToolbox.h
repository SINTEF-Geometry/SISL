/** ftFairingToolbox
//===========================================================================
//                                                                           
// File: ftFairingToolbox.h                                                   
//                                                                           
// Created: Wed Mar 07 2001                                         
//                                                                           
// Author: Vibeke Skytt 
//                                                                           
// Revision: 
//                                                                           
// Description: 
//
// Implementation files: ftFairingToolbox.C
// */                                                                          
//===========================================================================

#ifndef _FTFAIRINGTOOLBOX_H
#define _FTFAIRINGTOOLBOX_H

//===========================================================================
//===========================================================================

#include <vector>             // Standard library STL vector
#include <string>             // Standard library string
#include "tpTopologyTable.h"  // Topology container and analyzer
#include "ftEdgeBase.h"
#include "ftFaceBase.h"
#include "ftCurve.h"          // ftCurve returned by EdgeCurve(), Intersect()..
#include "ftPoint.h"
#include "ParamSurface.h"
#include "ftTangPriority.h"
#include "ftMessage.h"
#include "ftOffsetError.h"    // Contains eventual approximation errors
#include "ftSurfaceSet.h"
#include "ftSuperSurface.h"
#include "ftUtils.h"

// Forward declarations
//class ParamSurface;

typedef std::vector<std::shared_ptr<ParamSurface> > ftSurfaceGroup;
typedef std::vector<std::shared_ptr<SplineCurve> > ftCurveGroup;
typedef std::pair<double,double> surf_parameter;

/** A struct storing information about edge segments between
* surfaces. The ftBoundInt structs can be used to represent
* for instance gaps.
*/
struct ftBoundInt
{
  ParamSurface *surf1, *surf2;
  std::pair<surf_parameter, surf_parameter> parint1;
  std::pair<surf_parameter, surf_parameter> parint2;

  ftBoundInt() {surf1 = 0; surf2 = 0; }
  void addEntry(ParamSurface *s1, double par1_1[2], double par1_2[2],
		ParamSurface *s2, double par2_1[2], double par2_2[2])
  {
    surf1 = s1;
    surf2 = s2;
    parint1 = std::make_pair(std::make_pair(par1_1[0],par1_1[1]),
			     std::make_pair(par1_2[0],par1_2[1]));
    parint2 = std::make_pair(std::make_pair(par2_1[0],par2_1[1]),
			     std::make_pair(par2_2[0],par2_2[1]));
  }
};

/** A struct storing information about max and mean errors.
*/
struct ftApproxErrors
{
  double max_error_;
  double mean_error_;

  ftApproxErrors(double max, double mean)
    : max_error_(max), mean_error_(mean)
  {}
};

/** ftFairingToolbox -  Short description.
 * Detailed description.
 */
class ftFairingToolbox
{
public:
  
  // Constructors 
/** Constructor taking a subset of the topology tolerances as input. */
  ftFairingToolbox(double approxtol, double gap, double kink)
    : toptol_(gap, 10.0*gap, kink, 10.0*kink),
      approxtol_(approxtol),
      top_table_(gap, 10.0*gap, kink, 10.0*kink),
      sample_only_edges_(false) {}

  /** Constructor taking topology tolerances as input.  */
  ftFairingToolbox(double approxtol, double gap, double neighbour, 
		   double kink, double bend)
    : toptol_(tpTolerances(gap, neighbour, kink, bend)),
      approxtol_(approxtol),
      top_table_(gap, neighbour, kink, bend),
      sample_only_edges_(false) {}

  /// Destructor
  ~ftFairingToolbox();

    /** Set (or reset) topology tolerances The remaining tolerances
     *  are set depending on the given ones. */
    void setTolerances(double approxtol, double gap, double kink)
    {
	approxtol_ = approxtol;
	toptol_ = tpTolerances(gap, 10.0*gap, kink, 10.0*kink);
	top_table_.setTolerances(toptol_);
	
    }
    /** Set (or reset) topology tolerances */
    void setTolerances(double approxtol, double gap, double neighbour, 
		       double kink, double bend)
    {
	approxtol_ = approxtol;
	toptol_ = tpTolerances(gap, neighbour, kink, bend);
	top_table_.setTolerances(toptol_);
    }



  // Input functionality.
  /** Read geometry (surface, curves, bounded surfaces) in IGES format. */
  ftMessage readIges(std::istream& is);

  /** Read data points in offset format, approximate surfaces. */
  ftMessage readOffset(std::istream& is);

  /** Read surfaces from an application. Replaces any existing surfaces.*/
  ftMessage setSurfaces(std::vector<std::shared_ptr<ParamSurface> >& surfaces,
			std::vector<ftTangPriority>& surf_priority);

  /** Read groups of surfaces from an application. Replaces any existing surfaces.*/
  ftMessage setSurfaces(std::vector<ftSurfaceGroup>& surfaces,
			std::vector<ftTangPriority>& surf_priority);

  /** Read groups of lofting curves from an application. Replaces any existing surfaces.*/
  ftMessage setCurves(std::vector<ftCurveGroup >& curves,
		      std::vector<ftTangPriority>& surf_priority);

  /** Return constructed surfaces. */
  void getSurfaces(std::vector<std::shared_ptr<ParamSurface> >& surfaces);

  /** Return constructed surfaces and approximation errors. */
  void getSurfaces(std::vector<std::shared_ptr<ParamSurface> >& surfaces,
		   std::vector<ftApproxErrors>& approx_error);

  // Output functionality
  /** Write to IGES file */
  void writeIges(std::ostream& os);

    /** Obtain the approximation  errors from ReadOffset(). */
    const ftOffsetError& offsetError() { return offset_error_; }

   /** Construct the topology information regarding the input geometry */
  ftMessage buildTopology();

  /** Get information about holes. */
  int nmbHoles();

  /** The edges describing a hole is returned. */
  // 'whichhole' refers to element number in boundary_curves_,
  // not counting the outer boundaries.
  ftCurve getHole(int whichhole);

  /** Fill all hole(s) with surface(s). */
  ftMessage fillAllHoles(bool prefer_degenerate_surf);

  /** Fill a given hole with surfaces. */
  // 'whichhole' refers to element number in boundary_curves_,
  // not counting the outer boundaries.
  ftMessage fillHole(int whichhole, bool prefer_degenerate_surf = false);

  /** Check whether the hole given by ft_curve overlaps or not. */
  bool holeOverlaps(ftCurve& ft_curve);

  /** Return the number of pathwise disconnected objects. */
  int nmbObjects();

  /** Return the number of boundaries of a surface set (including holes). */
  int nmbBoundaries();

  /** Return a given boundary (may be a hole). */
  // 'whichhole' refers to element number in boundary_curves_.
  ftCurve getBoundary(int whichbound);

  /** Return information about all gaps. */
  ftCurve getGaps();

  /** Return information about all gaps. */
  void getGaps(std::vector<ftBoundInt>& gaps);

  /** Return information about all kinks. */
  ftCurve getKinks();

  /** Return information about all kinks. */
  void getKinks(std::vector<ftBoundInt>& kinks);

  /** Return information about all G1 discontinuities. */
  ftCurve getG1Disconts();

  /** Return information about all G1 discontinuities. */
  void getG1Disconts(std::vector<ftBoundInt>& g1disconts);

  // Manipulate the input geometry
  /** Update a surface set with regard to feature information, D1i. */
  ftMessage update();

  /** Remove gaps in the surface set */
  ftMessage removeGaps();

  /** Create surfaces from the given input data. */
  ftMessage createSurfaces();

/** Adapt expected continuity (mostly tangent plane) between neighbouring surfaces, D1i. */
  ftMessage adaptContinuity();

  /** Split all ftSurfaces in the fairing toolbox along G1 discontinuities. */
  void splitAlongKinks();

    /** Limiting the volume of interest. */
    void limitVolume(double xmin, double xmax,
		     double ymin, double ymax,
		     double zmin, double zmax);
    /** Check if a point is within the volume of interest. */
    bool pointWithinLimits(const ftPoint& point);

    /** Closest point. */
    ftPoint closestPoint(const ftPoint& point);
    /** Find all 'bad seams' or corners between neighbouring surfaces. */
    ftCurve cornerGapKinkCurves();
    /** Find all outer boundary curves. */
    ftCurve edgeCurve();
    /** Intersect the model with a plane. */
    ftCurve intersect(const ftPlane& plane);

    /** New functions to be implemented */
    /* Get information about the current limit volume */
    void getLimitVolume(double& xmin, double& xmax,
			double& ymin, double& ymax,
			double& zmin, double& zmax);

    /* Reset the limit volume to the complete model */
    void resetLimitVolume();

    /* Return information about T-connections within the box of interest */
    


  /** Connect input pts.
      The input pts contain ptr to par and spatial pt,
      as well as integer referring to input faces.
      If a 2d pts is missing the 3d pt is projected.
      If both d2 and 3d points exist we prefer 2d.
      If (interpolate == true) the 2d pts are interpolated,
      otherwise they are approximated.*/
  vector<ftCurve>
  connectPoints(std::vector<samplePoint>& pts,
		std::vector<std::shared_ptr<ftFaceBase> >& faces,
		double appr_tol,
		bool interpolate);


  /** Read groups of surfaces from an application. Replaces any existing surfaces.*/
  /** For a chart surface with less than 4 corners, user may define corners along a smooth edge.*/
  ftMessage setChartSurfaces(std::vector<ftSurfaceGroup>& surfaces,
			     std::vector<ftTangPriority>& surf_priority,
			     std::vector<std::pair<int, int> >& grid_res,
			     std::vector<std::vector<std::pair<Point, double> > >& edge_scales,
			     std::vector<bool>& symm_distr_functions,
			     std::vector<std::vector<Point> >& add_corner_pts);

  /** Make grid for all chart surfaces. */
  ftMessage makeGrid(RotationInfo* rot_info);

  /** Write grid for all chart surfaces, assuming sfs are made and
      grids generated. */
  ftMessage writeGrid(std::ofstream& os);

   void setFeatureTol(double featuretol);
  /** Set tolerance for approximation of feature curves. */

  void setSampleOnlyEdges()
    { sample_only_edges_ = true; }

  /** Read feature information from file, planned for D1i.*/
  ftMessage readFeature(std::istream& is);

  /** Read feature curves from an application, D1i. */
  ftMessage setFeature(std::vector<std::shared_ptr<SplineCurve> >& featurecurves);

   // -------------- Taken from ftGeometryModel --------------


  /** Read surfaces in Go format. */
  ftMessage readGo(std::istream& is);

   /** Constructs the actual topology. All members below should be called
     *  only after this member (it will work, but the results will be far
     *  less useful). */
    void constructTopology();

    /** Experimental, debug and development members. */
    Point boxmin() { return big_box_.low(); }
    Point boxmax() { return big_box_.high(); }
    void dumpSurfs(std::ostream& os);

protected:

// ------------- Some helper classes -------------------

class ftCell
{
protected:
    std::vector<ftSurface*> faces_;
    BoundingBox box_;
public:
    ftCell() {}
    ~ftCell() {}
    void setBox(const BoundingBox& box)
	{ box_ = box; }
    void addFace(ftSurface* f)
	{ faces_.push_back(f); }
    const BoundingBox& box() const
	{ return box_; }
    int num_faces() const
	{ return faces_.size(); }
    ftSurface* face(int i) const
        { return faces_[i]; }

};

class ftFaceInfo
{
public:
    ftSurface* face_;
    double dist_;
    bool closest_point_computed_;
    ftFaceInfo() : face_(0), dist_(-1), closest_point_computed_(false) {}
    ftFaceInfo(ftSurface* face, double dist, bool cpcomp) 
	: face_(face), dist_(dist), closest_point_computed_(cpcomp) {}
    bool operator < (const ftFaceInfo& fi) const { return dist_ < fi.dist_; }
};

class ftCellInfo
{
public:
    int index_;
    double dist_;
    ftCellInfo() :  index_(-1), dist_(-1) {}
    ftCellInfo(int i, double d) :  index_(i), dist_(d) {}
    bool operator < (const ftCellInfo& ci) const { return dist_ < ci.dist_; }
};


//------------------ Data members ----------------------
    tpTolerances toptol_;
    double approxtol_;
    double featuretol_;
    tpTopologyTable<ftEdgeBase, ftFaceBase> top_table_;
    std::vector<std::shared_ptr<ftFaceBase> > faces_;
    std::vector<ftApproxErrors> approx_errors_;
    bool sample_only_edges_;
    // For each separate object, we store first edge of all bnd curves.
    // First element is (what is supposed to be) the objects outer boundary.
    std::vector<std::vector<ftEdgeBase*> > boundary_curves_;
// Imported from ftGeometryModel
    std::vector<bool> face_checked_;
    int highest_face_checked_;
    BoundingBox big_box_;
    BoundingBox limit_box_;
    int ncellsx_, ncellsy_, ncellsz_;
    Point cell_delta_;
    std::vector<ftCell> cells_;
    std::vector<ftCellInfo> cell_info_;
    std::vector<ftPlane> limits_;
    ftOffsetError offset_error_;

    // ---------- Private functions ------------
 private:
    void addSegment(ftCurve& cv, ftEdgeBase* edge, ftCurveType ty);

    void getCurveofType(ftCurveType type, ftCurve& curve);
    void getCurveofType(ftCurveType type, 
			std::vector<ftBoundInt>& curve_intervals);
    // Set the elements in boundary_curves_ (based on top_table_).
    void setBoundaryCurves();

    // Get all edges meeting in corners in the inner of a surface set
    void getInnerCorners(std::vector<std::vector<ftEdgeBase*> >& edges_in_corner);

    // Get matching information along the current edge
    void getMatchingInfo(ftEdgeBase* startedge, std::vector<ftFaceBase* >& face,
			 std::vector<BoundaryPiece>& bdpiece);

// -------------- From ftGeometryModel ---------------------
    std::vector<ftCurveSegment> intersect(const ftPlane& plane, ftSurface* sf);
    ftCurve localIntersect(const ftPlane& plane, ftSurface* sf);
    void constructCells(int n1, int n2, int n3);
    void cellContaining(const Point& pt, int& i1, int& i2, int& i3);
    ftPoint closestPointLocal(const ftPoint& point);
    static double boxVecDist(const BoundingBox& b, const Point& v);
};

#endif // _FTFAIRINGTOOLBOX_H
      
       
  
