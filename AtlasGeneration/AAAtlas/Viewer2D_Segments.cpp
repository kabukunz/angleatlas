#include "Viewer2D_Segments.h"

#include <QKeyEvent>
#include <QWheelEvent>

#include "MeshViewer/OpenglHeaders.h"
#include <iostream>

Viewer2D_Segments::Viewer2D_Segments(QWidget *parent)
	: Viewer2D(parent)
{
	clear_widget();
}

Viewer2D_Segments::Viewer2D_Segments(QGLFormat& fmt, QWidget* parent /*= 0*/)
	: Viewer2D(fmt, parent)
{
	clear_widget();
}

void Viewer2D_Segments::clear_widget()
{
	segs.clear();
}

void Viewer2D_Segments::paintGL()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	setOrtho2D(width(), height());

	glColor3f(0.25, 0.25, 0.25);
	glBegin(GL_LINES);
	for (const auto& seg : mesh_edges)
	{
		glVertex2dv(seg.p0.data());
		glVertex2dv(seg.p1.data());
	}
	glEnd();

	glColor3f(1.0, 1.0, 1.0);
	glBegin(GL_LINES);
	for (const auto& seg : segs)
	{
		glVertex2dv(seg.p0.data());
		glVertex2dv(seg.p1.data());
	}
	glEnd();

// 	glColor3f(1.0, 0.0, 0.0);
// 	glPointSize(10.0f);
// 	glBegin(GL_POINTS);
// 	for (const auto& pt : pts)
// 	{
// 		glVertex2dv(pt.data());
// 	}
// 	glEnd();
}

void Viewer2D_Segments::adjust_range()
{
	if (segs.size() == 0)
	{
		return;
	}

	OpenMesh::Vec2d para_max, para_min;
	para_max = OpenMesh::Vec2d(-DBL_MAX, -DBL_MAX);
	para_min = OpenMesh::Vec2d(DBL_MAX, DBL_MAX);

	for (const auto& seg : segs)
	{
		para_max.maximize(seg.p0);
		para_min.minimize(seg.p0);

		para_max.maximize(seg.p1);
		para_min.minimize(seg.p1);
	}

	para_origin_0 = (para_max + para_min) * 0.5;
	para_range_0 = (para_max - para_min).max() * 0.55;

	para_origin = para_origin_0;
	para_range = para_range_0;
	updateGL();
}

void Viewer2D_Segments::add_segment(const OpenMesh::Vec2d& p0, const OpenMesh::Vec2d& p1)
{
	segs.emplace_back(segment{ p0, p1 });
}

void Viewer2D_Segments::add_point(const OpenMesh::Vec2d& p)
{
	pts.push_back(p);
}

void Viewer2D_Segments::add_point(double x, double y)
{
	pts.emplace_back(x, y);
}

void Viewer2D_Segments::add_mesh(const Mesh& m)
{
	for (auto e_h : m.edges())
	{
		auto p0 = m.point(m.to_vertex_handle(m.halfedge_handle(e_h, 0)));
		auto p1 = m.point(m.to_vertex_handle(m.halfedge_handle(e_h, 1)));
		mesh_edges.emplace_back(segment{ {p0[0], p0[1]}, {p1[0], p1[1]} });
	}
}
