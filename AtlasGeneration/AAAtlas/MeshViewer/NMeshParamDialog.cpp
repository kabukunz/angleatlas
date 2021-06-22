//#include "NMeshParamDialog.h"
#include <QApplication>
#include <QDesktopWidget>

#include "Common\CommonDefinitions.h" // HERE

#include "AAAtlas.h"
#include "MainViewerWidget.h"
#include <QSizePolicy>

MeshParamDialog::MeshParamDialog(QWidget* parent /* = 0 */)
	:QDialog(parent)
{
	initDialog();
}

MeshParamDialog::~MeshParamDialog()
{
}

void MeshParamDialog::initDialog()
{
	createWidget();
	createLayout();
}

void MeshParamDialog::createWidget()
{
	create_Basic_Operation_Information_Widget();
}

void MeshParamDialog::createLayout()
{
	tabWidget = new QTabWidget();

	QVBoxLayout *layout = new QVBoxLayout();
	layout->addWidget(tabWidget);
	setLayout(layout);

	setFixedWidth(220);
}

void MeshParamDialog::createTabs()
{
	create_info_tab();

	initTabs();
}

void MeshParamDialog::create_Basic_Operation_Information_Widget()
{
	QGridLayout* mainLayout = new QGridLayout();
	Basic_Operation_And_Information = new QWidget();
	Basic_Operation_And_Information->setLayout(mainLayout);

	view_BOI = new QScrollArea;
	view_BOI->setFocusPolicy(Qt::NoFocus);
	view_BOI->setFrameStyle(QFrame::NoFrame);
	view_BOI->setWidget(Basic_Operation_And_Information);
	view_BOI->setWidgetResizable(true);
}

void MeshParamDialog::create_info_tab()
{
	const auto& viewer = ((MainViewerWidget*)parent());
	const auto& MeshViewer = viewer->MeshViewer;
	const auto& mainWindow = (AAAtlas*)(viewer->parent());

	mesh_info_label = new QLabel();
	mesh_edge_label = new QLabel();
	mesh_face_label = new QLabel();

	QVBoxLayout *tabLayout = new QVBoxLayout();
	tabLayout->setAlignment(Qt::AlignTop);

	QVBoxLayout *infoLayout = new QVBoxLayout();
	infoLayout->setAlignment(Qt::AlignLeft | Qt::AlignTop);
	infoLayout->setSpacing(10);
	infoLayout->addWidget(mesh_info_label);
	infoLayout->addWidget(mesh_edge_label);
	infoLayout->addWidget(mesh_face_label);

	buttonPreprocess = new QPushButton(tr("Pre-Process"));
	buttonPolySquare = new QPushButton(tr("PolySquare"));
	buttonQuadCutting = new QPushButton(tr("Quad Cut"));
	buttonDisRuduce = new QPushButton(tr("Distortion Reduce"));

	QPushButton *buttonSigma = new QPushButton(tr("Sigma"));
	QPushButton *buttonPEBound = new QPushButton(tr("PE Bound"));

	labelSigma = new QLabel(tr("4.00"));
	labelPEBound = new QLabel(tr("0.8000"));
	labelGap = new QLabel(tr("On"));

	buttonGap = new QPushButton(tr("Gap"));
	buttonGap->setCheckable(true);



	connect(buttonPreprocess, &QPushButton::released, [&]() { MeshViewer->pre_process(); setPhase(1); });
	connect(buttonPolySquare, &QPushButton::released, [&]() { MeshViewer->polysquare(); setPhase(2); });
	connect(buttonQuadCutting, &QPushButton::released, [&]() { MeshViewer->quad_cutting(); setPhase(3); });
	connect(buttonDisRuduce, &QPushButton::released, [&]() { MeshViewer->distortion_reduce(); setPhase(4); });

	connect(buttonSigma, &QPushButton::released, [&]() { if (viewer->setParameters(tr("Kernel Width"), MeshViewer->kernel_width, 0.01, 1000.0, 2)) setSigmaLable(MeshViewer->kernel_width); });
	connect(buttonPEBound, &QPushButton::released, [&]() { if (viewer->setParameters(tr("PE Bound"), MeshViewer->pe_bound, 0.0, 0.99, 4)) setPEBound(MeshViewer->pe_bound); });
	connect(buttonGap, &QPushButton::released, [&]() { MeshViewer->chart_gap = (buttonGap->isChecked() ? 6.0 : 0.0); labelGap->setText(buttonGap->isChecked() ? tr("On") : tr("Off")); });

	QGridLayout *methodLayout = new QGridLayout();
	methodLayout->setAlignment(Qt::AlignTop);

	int iline = 0;
	methodLayout->addWidget(get_separator(), iline++, 0, 1, 2);

	methodLayout->addWidget(buttonSigma, iline, 0);
	methodLayout->addWidget(labelSigma, iline++, 1);
	methodLayout->addWidget(buttonPEBound, iline, 0);
	methodLayout->addWidget(labelPEBound, iline++, 1);
	methodLayout->addWidget(buttonGap, iline, 0);
	methodLayout->addWidget(labelGap, iline++, 1);

	methodLayout->addWidget(get_separator(), iline++, 0, 1, 2);

	methodLayout->addWidget(buttonPreprocess, iline++, 0, 1, 2);
	methodLayout->addWidget(buttonPolySquare, iline++, 0, 1, 2);
	methodLayout->addWidget(buttonQuadCutting, iline++, 0, 1, 2);
	methodLayout->addWidget(buttonDisRuduce, iline++, 0, 1, 2);

	tabLayout->addLayout(infoLayout);
	tabLayout->addLayout(methodLayout);
	tabWidget->setLayout(tabLayout);
}

QFrame* MeshParamDialog::get_separator()
{
	QFrame* line = new QFrame();
	line->setFrameShape(QFrame::HLine);
	line->setFrameShadow(QFrame::Sunken);
	return line;
}

void MeshParamDialog::initTabs()
{
	labelSigma->setText(tr("4.00"));
	labelPEBound->setText(tr("0.8000"));
	labelGap->setText(tr("On"));
	buttonGap->setChecked(true);

	updateMeshInfo();
	updateEdgeInfo(0.0, 0.0, 0.0);
	updateFaceInfo(0.0, 0.0, 0.0);

	setPhase(0);
	connect(((MainViewerWidget*)parent())->MeshViewer, &MeshViewerWidget::setPipelinePhase, this, &MeshParamDialog::setPhase);
}

void MeshParamDialog::updateMeshInfo()
{
	const auto& mesh = ((MainViewerWidget*)parent())->MeshViewer->mesh;

	QString mesh_info;
	mesh_info.sprintf("V: %d\nE: %d\nH: %d\nF: %d\nEuler: %d", mesh.n_vertices(), mesh.n_edges(), mesh.n_halfedges(), mesh.n_faces(), mesh.n_vertices() - mesh.n_edges() + mesh.n_faces());
	mesh_info_label->setText(mesh_info);
}

void MeshParamDialog::updateEdgeInfo(double min, double max, double avg)
{
	QString info;
	info.sprintf("Edge Length\nMax: %g\nMin: %g\nAvg: %g", max, min, avg);
	mesh_edge_label->setText(info);
}

void MeshParamDialog::updateFaceInfo(double min, double max, double avg)
{
	QString info;
	info.sprintf("Face Area\nMax: %g\nMin: %g\nAvg: %g", max, min, avg);
	mesh_face_label->setText(info);
}

void MeshParamDialog::setPhase(int i)
{
	buttonPreprocess->setEnabled(false);
	buttonPolySquare->setEnabled(false);
	buttonQuadCutting->setEnabled(false);
	buttonDisRuduce->setEnabled(false);

	switch (i)
	{
	case 0:
		buttonPreprocess->setEnabled(true);
		break;
	case 1:
		buttonPolySquare->setEnabled(true);
		break;
	case 2:
		buttonQuadCutting->setEnabled(true);
		break;
	case 3:
		buttonDisRuduce->setEnabled(true);
		break;
	}
}
