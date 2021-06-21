#ifndef VIEWER2D_PARAMETERIZATION_H
#define VIEWER2D_PARAMETERIZATION_H

#include <set>
#include <QGLWidget>
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include "MeshViewer/MeshDefinition.h"

#include "Viewer2D.h"

class Viewer2D_Parameterization : public Viewer2D
{
	Q_OBJECT

public:
	Viewer2D_Parameterization(const Mesh& _mesh, QWidget* parent = 0);
	Viewer2D_Parameterization(const Mesh& _mesh, QGLFormat& fmt, QWidget* parent = 0);

	virtual ~Viewer2D_Parameterization() = default;

	void clear_widget() override;
	void adjust_range() override;

	double get_bb_area()const { return bb_area; };

	std::set<int> selected_v, selected_e, selected_f;
	std::map<int, int> boundary_vk;

protected:
	const Mesh& mesh;

	void paintGL() override;
	void keyPressEvent(QKeyEvent *e);

	bool show_inner, show_k;
	double bb_area;
};

#endif // VIEWER2D_PARAMETERIZATION_H
