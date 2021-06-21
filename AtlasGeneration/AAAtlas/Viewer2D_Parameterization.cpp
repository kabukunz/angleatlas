#include "Viewer2D_Parameterization.h"

#include <QKeyEvent>
#include <QWheelEvent>

#include "MeshViewer/OpenglHeaders.h"
#include <iostream>

OpenMesh::Vec3uc color_pos[4] =
{
	{ 0,116,217 },		//1
	{ 57,204,204 },		//0
	{ 0, 250, 0 },		//-1
	{ 152, 56, 255 }	//-2-
};

OpenMesh::Vec3uc color_neg[4] =
{
	{ 255,65,54 },		//3
	{ 255,133,27 },		//4
	{ 255,220,0 },		//5
	{ 240,18,190 }		//6+
};

Viewer2D_Parameterization::Viewer2D_Parameterization(const Mesh& _mesh, QWidget *parent)
	: Viewer2D(parent), mesh(_mesh)
{
	clear_widget();
}

Viewer2D_Parameterization::Viewer2D_Parameterization(const Mesh& _mesh, QGLFormat& fmt, QWidget* parent /*= 0*/)
	: Viewer2D(fmt, parent), mesh(_mesh)
{
	clear_widget();
}

void Viewer2D_Parameterization::clear_widget()
{
	show_inner = true;
	show_k = false;
}

void Viewer2D_Parameterization::paintGL()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	setOrtho2D(width(), height());

	glColor3f(1.0, 0.5, 1.0);
	for (int f : selected_f)
	{
		glBegin(GL_POLYGON);
		for (auto fv : mesh.fv_range(mesh.face_handle(f)))
		{
			glVertex3dv(mesh.point(fv).data());
		}
		glEnd();
	}

	glBegin(GL_LINES);
	for (auto e_h : mesh.edges())
	{
		if (!mesh.is_boundary(e_h) && !show_inner) continue;

		glColor3ubv(mesh.color(e_h).data());

		auto h_h = mesh.halfedge_handle(e_h, 0);

		glVertex3dv(mesh.point(mesh.from_vertex_handle(h_h)).data());
		glVertex3dv(mesh.point(mesh.to_vertex_handle(h_h)).data());
	}
	glEnd();

	glPointSize(7.5f);
	if (show_k)
	{
		for (const auto& vk : boundary_vk)
		{
			int cur = vk.second;
			auto v_h = mesh.vertex_handle(vk.first);
			if (!mesh.is_boundary(v_h) || cur == 0) continue;

			const auto& color_list = (cur > 0) ? color_pos : color_neg;
			glColor3ubv(color_list[std::min(3, std::abs(cur) - 1)].data());
			glBegin(GL_POINTS);
			glVertex3dv(mesh.point(v_h).data());
			glEnd();
		}
	}

	glColor3f(1.0, 0.5, 0.0);
	glBegin(GL_POINTS);
	for (int v : selected_v)
	{
		glVertex3dv(mesh.point(mesh.vertex_handle(v)).data());
	}
	glEnd();

	glColor3f(1.0, 0.5, 1.0);
	glBegin(GL_LINES);
	for (int e : selected_e)
	{
		if (e == -1) continue;

		auto h_h = mesh.halfedge_handle(2 * e);

		glVertex3dv(mesh.point(mesh.from_vertex_handle(h_h)).data());
		glVertex3dv(mesh.point(mesh.to_vertex_handle(h_h)).data());
	}
	glEnd();
}

void Viewer2D_Parameterization::keyPressEvent(QKeyEvent *e)
{
	Viewer2D::keyPressEvent(e);

	switch (e->key())
	{
	case Qt::Key_I:
		show_inner = !show_inner;
		updateGL();
		break;
	case Qt::Key_K:
		show_k = !show_k;
		updateGL();
		break;
	}
}

void Viewer2D_Parameterization::adjust_range()
{
	if (mesh.vertices_empty())
	{
		return;
	}

	OpenMesh::Vec3d para_max, para_min;
	para_max = OpenMesh::Vec3d(-DBL_MAX, -DBL_MAX, -DBL_MAX);
	para_min = OpenMesh::Vec3d(DBL_MAX, DBL_MAX, DBL_MAX);

	for (auto v_h : mesh.vertices())
	{
		para_max.maximize(mesh.point(v_h));
		para_min.minimize(mesh.point(v_h));
	}

	OpenMesh::Vec3d origin3 = (para_max + para_min) * 0.5;
	para_origin_0[0] = origin3[0];
	para_origin_0[1] = origin3[1];

	auto para_diff = para_max - para_min;
	para_range_0 = para_diff.max() * 0.55;

	bb_area = para_diff[0] * para_diff[1];

	para_origin = para_origin_0;
	para_range = para_range_0;
	updateGL();
}
