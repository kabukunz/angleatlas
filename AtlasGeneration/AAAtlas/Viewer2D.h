#ifndef VIEWER2D_H
#define VIEWER2D_H

#include <QGLWidget>
#include <OpenMesh/Core/Geometry/VectorT.hh>

class Viewer2D : public QGLWidget
{
	Q_OBJECT

public:
	Viewer2D(QWidget* parent = 0);
	Viewer2D(QGLFormat& fmt, QWidget* parent = 0);

	virtual ~Viewer2D() = default;

	virtual void clear_widget() = 0;
	virtual void adjust_range() = 0;
	virtual void para_scale(double factor);

protected:
	virtual void initializeGL();
	virtual void paintGL() = 0;
	virtual void resizeGL(int width, int height);

	virtual void keyPressEvent(QKeyEvent *e);

	virtual void mousePressEvent(QMouseEvent *e);
	virtual void mouseMoveEvent(QMouseEvent *e);
	virtual void mouseReleaseEvent(QMouseEvent *e);
	virtual void wheelEvent(QWheelEvent * e);

	virtual void setOrtho2D(int w, int h);

	QPoint last_pos;
	OpenMesh::Vec2d para_origin, para_origin_0;
	double para_range, para_range_0;
};

#endif // VIEWER2D_H
