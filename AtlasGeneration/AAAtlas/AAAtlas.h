#pragma once

#include <QtWidgets/QMainWindow>
#include <QtGui>
#include <QtWidgets>

class QAction;
class MainViewerWidget;

class AAAtlas : public QMainWindow
{
	Q_OBJECT

public:
	AAAtlas(QWidget *parent = Q_NULLPTR);
	~AAAtlas();

	void open_mesh_from_main(char* filemane);

	MainViewerWidget* viewer;

protected:

private slots:
	bool save();
	bool saveAs();

	void setMouseMode_slot(int mm);
	void setDrawMode_slot(int dm);

private:
	void initWindow();
	void createActions();
	void createMenus();
	void createToolBars();
	void createStatusBar();

private slots:
	void wireFrameShow();
	void hiddenLinesShow();
	void solidFlatShow();
	void flatPointsShow();
	void solidSmoothShow();
	void pointSetShow();
	void curvatureShow();
	void DrawBBox();
	void DrawMeshBoundary();
	void pointPick();
	void vertexPick();
	void edgePick();
	void facePick();
	void moveVertex();
	void edit_undo();
	void set_edit_undo_enable(bool b);
	void edit_redo();
	void set_edit_redo_enable(bool b);
	void edge_collpase();
	void edge_flip();
	void edge_split();

	void aux_inverse_mesh_connectivity();
	void aux_scale_mesh_BBox();
	void aux_split_quad_mesh();
	void aux_transform_mesh();
	void aux_generate_rect_mesh();

	void aux_find_vertex_by_index();
	void aux_find_edge_by_index();
	void aux_find_face_by_index();
	void aux_find_vertex_by_valance();

	void aux_delete_vertex_valence_four();
	void aux_delete_vertex_valence_three();
	void aux_split_vertex_valence_eight();

	//my
	void aux_my_process();
	//

private:
	QAction* switch_tools;

	// File Actions.
	QAction* openAction;
	QAction* saveAction;
	QAction* saveParaAction;
	QAction* saveAsAction;
	QAction* exitAction;
	QAction* clearAction;

	// Edit Actions.
	QAction* updateMeshAction;

	// View Actions.
	QAction* wireFrameAction;
	QAction* hiddenLinesAction;
	QAction* solidFlatAction;
	QAction* flatPointsAction;
	QAction* solidSmoothAction;
	QAction* pointSetAction;
	QAction* curvatureAction;

	// other View Actions.
	QAction* drawBoundingBox;
	QAction* drawMeshBoundaryAction;

	//auxiliary operation
	QAction* inverse_mesh_connectivity;
	QAction* scale_mesh_BBox;
	QAction* split_quad_mesh;
	QAction* transform_mesh;
	QAction* generate_rect_mesh;
	QAction* find_face_by_id;
	QAction* find_edge_by_id;
	QAction* find_vertex_by_id;
	QAction* find_vertex_by_valance;
	QAction* delete_vertex_valence_four;
	QAction* delete_vertex_valence_three;
	QAction* split_vertex_valence_eight;

	//my
	QAction* myprocess;

	//mouse Actions
	QAction* PointPickAction;
	QAction* VertexPickAction;
	QAction* FacePickAction;
	QAction* EdgePickAction;
	QAction* clearSelectedAction;
	QAction* moveVertexAction;

	QAction* EditUndoAction;
	QAction* EditRedoAction;
	QAction* EdgeCollpaseAction;
	QAction* EdgeFlipAction;
	QAction* EdgeSplitAction;

	// Help Actions.
	QAction* aboutAction;
	QAction* aboutQtAction;

	QAction* squareWidgetAction;
	QAction* saveScreenAction;
	QAction* saveScreenActionDirectly;

	// Menus.
	QMenu* fileMenu;
	QMenu* viewMenu;
	QMenu* mouseMenu;
	QMenu* Auxiliary_Menu;
	QMenu* helpMenu;
	QMenu* toolsMenu;

	// ToolBars.
	QToolBar* fileToolBar;
	QToolBar* viewToolBar;
	QToolBar* otherViewBar;
	QToolBar* mouseToolBar;
	QToolBar* localOperationBar;
	QToolBar* findByIdBar;
	QToolBar* saveToolBar;

	//label
	QLabel* statusLabel;

private:
	void setAllMouseActionChecked(bool b);
	void setAllViewActionChecked(bool b);

private:
	int latestDrawMode;

};
