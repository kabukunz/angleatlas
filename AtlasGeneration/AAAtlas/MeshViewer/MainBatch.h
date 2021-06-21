#ifndef MESHPROCESSING_MAIN_VIEWWE_WIDGET_H
#define MESHPROCESSING_MAIN_VIEWWE_WIDGET_H

#include <QtGui>
#include <QString>
#include <QMessageBox>
#include <QFileDialog>
//main widget
#include "InteractiveViewerWidget.h"
#include "MeshParamDialog.h"

class MainViewerWidget : public QSplitter
{
	Q_OBJECT
public:
	MainViewerWidget(QWidget* _parent=0);
	~MainViewerWidget();

	InteractiveViewerWidget* MeshViewer;

	void setDrawMode(int dm)
	{
		MeshViewer->setDrawMode(dm);
	}
	void setMouseMode(int dm)
	{
		MeshViewer->setMouseMode(dm);
	}

	void set_show_BBox()
	{
		MeshViewer->set_draw_bbox_ok();
	}
	void set_show_mesh_boundary()
	{
		MeshViewer->set_draw_mesh_boundary_ok();
	}
	void openMesh_fromMain(char* filename)
	{
		std::cout << "Open Mesh" << std::endl;
		QString str(filename);
		open_mesh_gui(filename);
	}

	void edit_undo()
	{
		MeshViewer->edit_undo_viewer();
	}

	void edit_redo()
	{
		MeshViewer->edit_redo_viewer();
	}

	void setParamVisible(bool t)
	{
		MeshParam->setVisible(t);
	}

	bool setParameters(const QString& title, double& p, double min, double max, int decimals);

public slots:

	void open_mesh_query()
	{
		QString fileName = QFileDialog::getOpenFileName(this,
			tr("Open mesh file"),
			lastpath,
			tr("OBJ Files (*.obj);;"
			"OFF Files (*.off);;"
			"PLY Files (*.ply);;"
			"STL Files (*.stl);;"
			"All Files (*)"));
		if (!fileName.isEmpty())
		{
			open_mesh_gui(fileName);
		}
	}
	void save_mesh_query() 
	{
		QString fileName = QFileDialog::getSaveFileName(this,
			tr("Save mesh file"),
			lastpath,
			tr("OBJ Files (*.obj);;"
			"OFF Files (*.off);;"
			"PLY Files (*.ply);;"
			"STL Files (*.stl);;"
			"All Files (*)"));
		if (!fileName.isEmpty())
		{
			save_mesh_gui(fileName);
		}
	}
	void saveOpenGLScreen()
	{
		QString fileName = QFileDialog::getSaveFileName(this,
			("Save screen as image file"),
			("../Results/untitled.png"),
			("PNG Files (*.png);;BMP Files (*.bmp);;JPG Files (*.jpg);;"
			"All Files (*)"));
		if (!fileName.isEmpty())
		{
			save_screen_gui(fileName);
		}
	}
	void save_opengl_screen(const QString& str)
	{
		MeshViewer->saveScreen(str.toLocal8Bit());
	}
	virtual void update_mesh()
	{
		if( MeshViewer->mesh_ref().n_vertices() != 0 )
		{
			MeshViewer->updateMesh();
		}
	}

	virtual void clear_all_mesh()
	{
		if(LoadMeshSuccess)
		{
			LoadMeshSuccess = false;
			MeshViewer->clearAllMesh();
			MeshParam->initTabs();
		}
	}

	virtual void clear_all_selected()
	{
		if(LoadMeshSuccess)
		{
			MeshViewer->clearSelectedData();
			MeshViewer->updateGL();
		}
	}

	void LoadMeshFromInner(bool OK, const QString& fname)
	{
		LoadMeshSuccess = OK;
		if(LoadMeshSuccess)
		{
			SetMeshForALL();
			lastpath = QFileInfo(fname).absolutePath();
		}
		emit( haveLoadMesh(fname) );
	};

signals:
	void haveLoadMesh(const QString& filePath);
	void setMouseMode_signal_main(int);
	void setDrawMode_signal_main(int);

	void set_edit_undo_enable_signal(bool);
	void set_edit_redo_enable_signal(bool);

protected:
	virtual void initViewerWindow();
	virtual void createParamIDialog();
	virtual void createViewerDialog();
	virtual void save_mesh_gui(const QString& fname);
	virtual void open_mesh_gui(const QString& fname);
	virtual void save_screen_gui(const QString& fname);
	virtual void load_file(bool success, const QString& fileName);

protected:
	bool LoadMeshSuccess;

	QString lastpath;

public:
	MeshParamDialog * MeshParam;

private:
	
	void SetMeshForALL( )
	{
	}

#pragma region Auxiliary Function
public:
	void aux_inverse_mesh_connectivity()
	{
		MeshViewer->inverse_mesh_connectivity();
	}

	void aux_scale_mesh_using_BBox(int max_len)
	{
		MeshViewer->scale_mesh_using_BBox(max_len);
	}

	void aux_split_quad_mesh()
	{
		MeshViewer->split_quad_mesh();
	}

	void transform_mesh(const std::vector<double>& m)
	{
		MeshViewer->transform_mesh(m);
	}

	void generate_rect_mesh(const std::vector<double>& m)
	{
		MeshViewer->geneate_rect_mesh(m);
	}

	void aux_find_vertex_by_id(int id)
	{
		MeshViewer->find_vertex_by_id(id);
	}

	void aux_find_face_by_id(int id)
	{
		MeshViewer->find_face_by_id(id);
	}

	void aux_find_edge_by_id(int id)
	{
		MeshViewer->find_edge_by_id(id);
	}

	void aux_find_vertex_by_valance(int valance)
	{
		MeshViewer->find_vertex_by_valance( valance );
	}

	void aux_delete_vertex_valence_four()
	{
		MeshViewer->delete_vertex_valence_four( );
	}

	void aux_delete_vertex_valence_three()
	{
		MeshViewer->delete_vertex_valence_three( );
	}

	void aux_split_vertex_valence_eight()
	{
		MeshViewer->split_vertex_valence_eight( );
	}

	//my
	void aux_my_process()
	{
		std::vector<std::string> files;         //保存该目录下所有文件的文件名
		__int64 hFile = 0;
		struct __finddata64_t fileinfo;
		hFile = _findfirst64(".\\input6\\*.obj", &fileinfo);
		while (!_findnext64(hFile, &fileinfo))
		{
			files.push_back(fileinfo.name);
		}
		_findclose(hFile);

		for (int i = 1; i < files.size(); i++)
		{
			QString in_file_name = QString::fromStdString(".\\input6\\"+files[i]);
			open_mesh_gui(in_file_name.toLocal8Bit());
			std::cout << files[i] << std::endl;
			/*Mesh input_mesh = MeshViewer->mesh;
			bool my_is_flip = 0;
			for (auto f_h : input_mesh.faces())
			{
				std::vector<OpenMesh::Vec3d> face_point;
				for (auto fv_h : input_mesh.fv_range(f_h))
				{
					face_point.push_back(input_mesh.point(fv_h));
				}
				OpenMesh::Vec3d vec1 = face_point[1] - face_point[0];
				OpenMesh::Vec3d vec2 = face_point[2] - face_point[0];
				double my_det_p = OpenMesh::cross(vec1, vec2).norm();
				if (my_det_p <= 0)
				{
					my_is_flip = 1;
					break;
				}
			}
			if (my_is_flip)
			{
				continue;
			}*/

			std::string out_file_name = files[i].substr(0, files[i].find_last_of('.')) + "_axis.obj";
			MeshViewer->my_file_name = out_file_name;

			MeshViewer->M_PE = 0.8;
			MeshViewer->pre_process();
			MeshViewer->polysquare();
			MeshViewer->quad_cutting();
			MeshViewer->distortion_reduce();
		}
	}
#pragma endregion

};


#endif