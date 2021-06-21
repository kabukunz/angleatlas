#include "MainViewerWidget.h"
#include <QDir>
#include <fstream>
#include <streambuf>
#include <stdio.h>
#include <queue>

#include "Mesh_doubleIO.h" // HERE

MainViewerWidget::MainViewerWidget(QWidget* _parent/* =0 */)
{
	initViewerWindow();
	LoadMeshSuccess = false;

	lastpath = tr(".");
}
MainViewerWidget::~MainViewerWidget()
{
};

void MainViewerWidget::initViewerWindow()
{
	createParamIDialog();
	createViewerDialog();

	this->setOrientation( Qt::Horizontal );

	//this->addWidget(debugDialog);
	//OpenGL mesh viewer
	this->addWidget(MeshParam);
	this->addWidget(MeshViewer);

	//set the splitter line color
	//this->setStyleSheet("QSplitter::handle { background-color: green }");
	//QSplitterHandle* splitterHandle = this->handle(1);
	//splitterHandle->setDisabled(true);

	connect(MeshViewer,SIGNAL(setMouseMode_signal(int)),SIGNAL(setMouseMode_signal_main(int)));
	connect(MeshViewer,SIGNAL(setDrawMode_signal(int)),SIGNAL(setDrawMode_signal_main(int)));
	connect(MeshViewer,SIGNAL(set_edit_undo_enable_viewer_signal(bool)),SIGNAL(set_edit_undo_enable_signal(bool)));
	connect(MeshViewer,SIGNAL(set_edit_redo_enable_viewer_signal(bool)),SIGNAL(set_edit_redo_enable_signal(bool)));
}

void MainViewerWidget::createParamIDialog()
{
	MeshParam = new MeshParamDialog();
}

void MainViewerWidget::createViewerDialog()
{
	QGLFormat glFormat;
	glFormat.setSampleBuffers(true);
	glFormat.setSamples(16);

	MeshViewer = new InteractiveViewerWidget(glFormat, NULL);
	MeshViewer->setAcceptDrops(true);
	connect(MeshViewer,SIGNAL(loadMeshOK(bool,QString)), this, SLOT(LoadMeshFromInner(bool,QString)) );
	//connect(MeshViewer,SIGNAL(clear_your_data()), this, SLOT(clear_all_data()) );
}

void MainViewerWidget::open_mesh_gui(const QString& fname)
{
	if (fname.isEmpty() || !MeshViewer->openMesh(fname.toLocal8Bit())) 
	{
		QString msg = "Cannot read mesh from file:\n '";
		msg += fname;
		msg += "'";
		QMessageBox::critical(NULL, windowTitle(), msg);
	}
	else
	{
		LoadMeshSuccess = true;
		lastpath = QFileInfo(fname).absolutePath();
		MeshViewer->setDrawMode(InteractiveViewerWidget::FLAT_POINTS);
		MeshViewer->setMouseMode(InteractiveViewerWidget::TRANS);
		if(LoadMeshSuccess)
		{
			SetMeshForALL();
		}
		emit(haveLoadMesh(fname));
		MeshViewer->load_parameterization();
		MeshViewer->input_distortion = MeshViewer->para_cutting->get_distortion();
		MeshViewer->updateGL();
	}
}

void MainViewerWidget::save_mesh_gui(const QString& fname)
{
	if (fname.isEmpty() || !MeshViewer->saveMesh(fname.toLocal8Bit()))
	{
		QString msg = "Cannot read mesh from file:\n '";
		msg += fname;
		msg += "'";
		QMessageBox::critical(NULL, windowTitle(), msg);
	}
}

void MainViewerWidget::save_screen_gui(const QString& fname)
{
	if (fname.isEmpty() || !MeshViewer->saveScreen(fname.toLocal8Bit()))
	{
		QString msg = "Cannot save image to file:\n '";
		msg += fname;
		msg += "'";
		QMessageBox::critical(NULL, windowTitle(), msg);
	}
}

void MainViewerWidget::load_file(bool success, const QString& fileName)
{
	if (!success) 
	{
		QMessageBox::critical(NULL, windowTitle(), tr("Cannot read file."));
	}
	else
	{
		MeshViewer->setDrawMode(InteractiveViewerWidget::FLAT_POINTS);
		MeshViewer->setMouseMode(InteractiveViewerWidget::TRANS);

		lastpath = QFileInfo(fileName).absolutePath();
	}
}

bool MainViewerWidget::setParameters(const QString& title, double& p, double min, double max, int decimals)
{
	bool isOK;
	double new_p = QInputDialog::getDouble(nullptr, title, title, p, min, max, decimals, &isOK, Qt::WindowFlags() | Qt::MSWindowsFixedSizeDialogHint);

	if (isOK) p = new_p;
	return isOK;
}