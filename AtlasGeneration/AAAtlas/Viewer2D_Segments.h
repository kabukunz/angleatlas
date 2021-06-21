#ifndef VIEWER2D_SEGMENTS_H
#define VIEWER2D_SEGMENTS_H

#include <QGLWidget>
#include <OpenMesh/Core/Geometry/VectorT.hh>

#include "Viewer2D.h"
#include "MeshViewer/MeshDefinition.h"

class Viewer2D_Segments : public Viewer2D
{
	Q_OBJECT

public:
	Viewer2D_Segments(QWidget* parent = 0);
	Viewer2D_Segments(QGLFormat& fmt, QWidget* parent = 0);

	virtual ~Viewer2D_Segments() = default;

	void clear_widget() override;
	void adjust_range() override;

	void add_segment(const OpenMesh::Vec2d& p0, const OpenMesh::Vec2d& p1);
	void add_point(const OpenMesh::Vec2d& p);
	void add_point(double x, double y);
	void add_mesh(const Mesh& m);

protected:
	void paintGL() override;
	
	struct segment
	{
		OpenMesh::Vec2d p0, p1;
	};

	std::vector<segment> segs;
	std::vector<segment> mesh_edges;
	std::vector<OpenMesh::Vec2d> pts;
};

#endif // VIEWER2D_SEGMENTS_H
