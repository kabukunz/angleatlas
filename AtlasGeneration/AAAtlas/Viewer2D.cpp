#include "Viewer2D.h"

#include <QKeyEvent>
#include <QWheelEvent>

#include "MeshViewer/OpenglHeaders.h"
#include <iostream>

Viewer2D::Viewer2D(QWidget *parent)
	: QGLWidget(parent)
{
}

Viewer2D::Viewer2D(QGLFormat& fmt, QWidget* parent /*= 0*/)
	: QGLWidget(fmt, parent)
{
}

void Viewer2D::para_scale(double factor)
{
	para_origin *= factor;
	para_origin_0 *= factor;
	para_range *= factor;
	para_range_0 *= factor;
}

void Viewer2D::initializeGL()
{
	setGeometry(200, 200, 800, 800);
	glShadeModel(GL_FLAT);
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClearDepth(1.0);

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
}

void Viewer2D::resizeGL(int width, int height)
{
	glLoadIdentity();

	glViewport(0, 0, width, height);
	setOrtho2D(width, height);
}

void Viewer2D::keyPressEvent(QKeyEvent *e)
{
	switch (e->key())
	{
	case Qt::Key_R:
		para_origin = para_origin_0;
		para_range = para_range_0;
		updateGL();
		break;
	case Qt::Key_Escape:
		close();
	}
}

void Viewer2D::mousePressEvent(QMouseEvent *e)
{
	if (e->button() == Qt::LeftButton)
	{
		last_pos = e->pos();
		setCursor(Qt::SizeAllCursor);
	}
}

void Viewer2D::mouseMoveEvent(QMouseEvent *e)
{
	if (e->buttons() & Qt::LeftButton)
	{
		QPoint win_move = e->pos() - last_pos;
		int win_range = std::min(width(), height());
		OpenMesh::Vec2d para_move(win_move.x(), -win_move.y());
		para_move *= 2 * para_range / win_range;
		para_origin -= para_move;

		last_pos = e->pos();
		updateGL();
	}
}

void Viewer2D::mouseReleaseEvent(QMouseEvent *e)
{
	if (e->button() == Qt::LeftButton)
	{
		last_pos = e->pos();
		setCursor(Qt::ArrowCursor);
	}
}

void Viewer2D::wheelEvent(QWheelEvent * e)
{
	para_range *= std::pow(0.999, e->delta());
	updateGL();
}

void Viewer2D::setOrtho2D(int w, int h)
{
	double range_x, range_y;

	if (w < h)
	{
		range_x = para_range;
		range_y = (double)h / w * para_range;
	}
	else
	{
		range_x = (double)w / h * para_range;
		range_y = para_range;
	}

	gluOrtho2D(para_origin[0] - range_x, para_origin[0] + range_x, para_origin[1] - range_y, para_origin[1] + range_y);
}