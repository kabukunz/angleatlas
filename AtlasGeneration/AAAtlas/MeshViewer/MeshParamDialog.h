#ifndef MESHPROCESSING_MESHPARAMDIALOG_H
#define MESHPROCESSING_MESHPARAMDIALOG_H

#include <QDialog>
#include <QtGui>
#include <QtWidgets>
#include <iostream>
#include <string>

class MeshParamDialog : public QDialog
{
	Q_OBJECT
public:
	MeshParamDialog(QWidget* parent = 0);
	~MeshParamDialog();

	QSize sizeHint()
	{
		QRect rect = QApplication::desktop()->screenGeometry();
		return QSize(int(rect.width()*0.15), rect.height());
	}

	void createTabs();
	void initTabs();

	void updateMeshInfo();
	void updateEdgeInfo(double min, double max, double avg);
	void updateFaceInfo(double min, double max, double avg);

	void setSigmaLable(double x)
	{
		labelSigma->setText(QString::number(x, 'f', 2));
	}
	void setPEBound(double x)
	{
		labelPEBound->setText(QString::number(x, 'f', 4));
	}

	void setPhase(int i);

private:
	QTabWidget * tabWidget;

	QWidget* tabInfo;
	QLabel *mesh_info_label;
	QLabel *mesh_edge_label;
	QLabel *mesh_face_label;

	QLabel *labelSigma;
	QLabel *labelPEBound;
	QLabel *labelGap;
	QPushButton *buttonGap = nullptr;

	QWidget* tabMod;

	QPushButton *buttonPreprocess = nullptr;
	QPushButton *buttonPolySquare = nullptr;
	QPushButton *buttonQuadCutting = nullptr;
	QPushButton *buttonDisRuduce = nullptr;

signals:

private:
	QWidget * Basic_Operation_And_Information;
	QScrollArea *view_BOI;

	QLabel* leftLabel_BOI;

private:
	void create_Basic_Operation_Information_Widget();

private:
	void initDialog();
	void createWidget();
	void createLayout();

	void create_info_tab();

	QFrame* get_separator();
};

#endif