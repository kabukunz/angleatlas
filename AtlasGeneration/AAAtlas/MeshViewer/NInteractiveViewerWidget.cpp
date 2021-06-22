// #include <QMouseEvent>
// #include <QLineEdit>
// #include <QDragEnterEvent>
// #include <QDropEvent>
// #include <QtCore>
// #include <QUrl>
// #include <QDir>
// #include <QInputDialog>
#include <queue>

#include "Mesh_doubleIO.h" // HERE

#include "Eigen/Eigen" // HERE

#include "NInteractiveViewerWidget.h"
#include "Scaffold/StateManager.h" // HERE

//InteractiveViewerWidget::InteractiveViewerWidget(QWidget* parent /* = 0 */)
//	:MeshViewerWidget(parent)


NInteractiveViewerWidget::NInteractiveViewerWidget() : NMeshViewerWidget ()
{
	draw_new_mesh = false;
	// clearSelectedData();
	kdTree = NULL;

	kernel_width = 4.0;
	goal_factor = 0.1;
	energy_exp_factor = 1.0;
	split_thres = 0.02;
}

// InteractiveViewerWidget::InteractiveViewerWidget(QGLFormat& _fmt, QWidget* _parent)
// :MeshViewerWidget(_fmt, _parent)
// {
// 	draw_new_mesh = false;
// 	clearSelectedData();
// 	kdTree = NULL;

// 	kernel_width = 4.0;
// 	goal_factor = 0.1;
// 	energy_exp_factor = 1.0;
// 	split_thres = 0.02;
// }

NInteractiveViewerWidget::~NInteractiveViewerWidget()
{
	if(kdTree) delete kdTree;
}

void NInteractiveViewerWidget::initMesh()
{
	NMeshViewerWidget::initMesh();

	// selectedVertex = new_boundary_vertices;
}

// void InteractiveViewerWidget::setMouseMode(int mm)
// {
// 	if(mouse_mode_ != T2_MODE)
// 	{
// 		mouse_mode_ = mm;
// 		if( TRANS != mouse_mode_ )
// 		{ buildIndex(); }
// 		emit setMouseMode_signal(mm);
// 	}
// }

// void InteractiveViewerWidget::mousePressEvent(QMouseEvent *_event)
// {
// 	if(mouse_mode_ == TRANS)
// 	{
// 		MeshViewerWidget::mousePressEvent(_event);
// 	}
// 	else
// 	{
// 		if(mouse_mode_ != T2_MODE)
// 		{
// 			pick_point( _event->x(), _event->y() );
// 			if(mouse_mode_ == VERTEXPICK)
// 			{
// 				pick_vertex( _event->x(), _event->y() );
// 			}
// 			else if(mouse_mode_ == FACEPICK)
// 			{
// 				pick_face( _event->x(), _event->y() );
// 			}
// 			else if(mouse_mode_ == EDGEPICK)
// 			{
// 				pick_edge( _event->x(), _event->y() );
// 			}
// 			else if(mouse_mode_ == POINTPICK)
// 			{
// 			}
// 			else if( mouse_mode_ == MOVE )
// 			{
// 				pick_vertex( _event->x(), _event->y() );//set the selected handle
// 			}
// 			else if(mouse_mode_ == EDGECOLLAPSE)
// 			{
// 				int desired_edge = find_edge_using_selected_point();
// 				if(desired_edge >= 0) 
// 				{
// 					Mesh::HalfedgeHandle heh = mesh.halfedge_handle( mesh.edge_handle(desired_edge), 0 );
// 					OpenMesh::Vec3d from_p = mesh.point(mesh.from_vertex_handle(heh));
// 					OpenMesh::Vec3d to_p = mesh.point(mesh.to_vertex_handle(heh));
// 					OpenMesh::Vec3d sp(selectedPoint[0],selectedPoint[1],selectedPoint[2]);
// 					bool collapse_ok = true;
// 					if( (sp-from_p).sqrnorm() > (to_p-sp).sqrnorm() )
// 					{
// 						if( mesh.is_collapse_ok(heh) )
// 						{
// 							mesh.collapse(heh);
// 						}
// 						else
// 						{
// 							collapse_ok = false;
// 							printf("[%d] Collapse Not OK!\n", desired_edge);
// 						}
// 					}
// 					else
// 					{
// 						heh = mesh.opposite_halfedge_handle(heh);
// 						if( mesh.is_collapse_ok(heh) )
// 						{
// 							mesh.collapse(heh);
// 						}
// 						else
// 						{
// 							collapse_ok = false;
// 							printf("[%d] Collapse Not OK!\n", desired_edge);
// 						}
// 					}
// 					if(collapse_ok)
// 					{
// 						mesh.garbage_collection();
// 						buildIndex();
// 						if( mesh_vector.size() - 1 > mesh_vector_index )
// 						{
// 							mesh_vector.erase( mesh_vector.begin() + mesh_vector_index + 1, mesh_vector.end() );
// 						}
// 						mesh_vector.push_back( mesh ); mesh_vector_index += 1;
// 						emit set_edit_undo_enable_viewer_signal( true );
// 						emit set_edit_redo_enable_viewer_signal( false );
// 					}
// 					clearSelectedData();
// 				}
// 			}
// 			else if (mouse_mode_ == EDGEFLIP)
// 			{
// 				int desired_edge = find_edge_using_selected_point();
// 				if(desired_edge >= 0) 
// 				{
// 					Mesh::EdgeHandle eh = mesh.edge_handle(desired_edge);
// 					if( is_flip_ok_openmesh(eh, mesh))
// 					{
// 						flip_openmesh(eh, mesh);
// 						if( mesh_vector.size() - 1 > mesh_vector_index )
// 						{
// 							mesh_vector.erase( mesh_vector.begin() + mesh_vector_index + 1, mesh_vector.end() );
// 						}
// 						mesh_vector.push_back( mesh ); mesh_vector_index += 1;
// 						emit set_edit_undo_enable_viewer_signal( true );
// 						emit set_edit_redo_enable_viewer_signal( false );
// 					}
// 					else
// 					{
// 						printf("[%d] Flip Not OK!\n", desired_edge);
// 					}
// 					clearSelectedData();
// 				}
// 			}
// 			else if (mouse_mode_ == EDGESPLIT)
// 			{
// 				int desired_edge = find_edge_using_selected_point();
// 				if(desired_edge >= 0) 
// 				{
// 					Mesh::EdgeHandle eh = mesh.edge_handle(desired_edge);
// 					Mesh::HalfedgeHandle heh = mesh.halfedge_handle( eh, 0 );
// 					Mesh::HalfedgeHandle heh_ = mesh.halfedge_handle( eh, 1 );
// 					Mesh::VertexHandle vh0 = mesh.to_vertex_handle(heh);
// 					Mesh::VertexHandle vh1 = mesh.to_vertex_handle(heh_);
// 					OpenMesh::Vec3d s = mesh.point( vh1 );
// 					OpenMesh::Vec3d e = mesh.point( vh0 );
// 					Mesh::VertexHandle vh = mesh.add_vertex( (s + e)*0.5 );
// 					std::vector<Mesh::VertexHandle> one_face(3);
// 					if(mesh.is_boundary(eh))
// 					{
// 						if(Mesh::InvalidFaceHandle != mesh.face_handle(heh))
// 						{
// 							Mesh::VertexHandle vh2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(heh));
// 							mesh.delete_edge(eh, false); mesh.garbage_collection();
// 							one_face[0] = vh0; one_face[1] = vh2; one_face[2] = vh; mesh.add_face(one_face);
// 							one_face[0] = vh2; one_face[1] = vh1; one_face[2] = vh; mesh.add_face(one_face);
// 						}
// 						else
// 						{
// 							Mesh::VertexHandle vh3 = mesh.to_vertex_handle(mesh.next_halfedge_handle(heh_));
// 							mesh.delete_edge(eh, false); mesh.garbage_collection();
// 							one_face[0] = vh3; one_face[1] = vh0; one_face[2] = vh; mesh.add_face(one_face);
// 							one_face[0] = vh1; one_face[1] = vh3; one_face[2] = vh; mesh.add_face(one_face);
// 						}
// 					}
// 					else
// 					{
// 						Mesh::VertexHandle vh2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(heh));
// 						Mesh::VertexHandle vh3 = mesh.to_vertex_handle(mesh.next_halfedge_handle(heh_));
// 						mesh.delete_edge(eh, false); mesh.garbage_collection();
// 						one_face[0] = vh0; one_face[1] = vh2; one_face[2] = vh; mesh.add_face(one_face);
// 						one_face[0] = vh2; one_face[1] = vh1; one_face[2] = vh; mesh.add_face(one_face);
// 						one_face[0] = vh3; one_face[1] = vh0; one_face[2] = vh; mesh.add_face(one_face);
// 						one_face[0] = vh1; one_face[1] = vh3; one_face[2] = vh; mesh.add_face(one_face);
// 					}

// 					mesh.update_normals();
// 					buildIndex();
// 					clearSelectedData();

// 					if( mesh_vector.size() - 1 > mesh_vector_index )
// 					{
// 						mesh_vector.erase( mesh_vector.begin() + mesh_vector_index + 1, mesh_vector.end() );
// 					}
// 					mesh_vector.push_back( mesh ); mesh_vector_index += 1;
// 					emit set_edit_undo_enable_viewer_signal( true );
// 					emit set_edit_redo_enable_viewer_signal( false );
// 				}
// 			}
// 		}
// 	}
// 	updateGL();
// }

// void InteractiveViewerWidget::mouseMoveEvent(QMouseEvent *_event)
// {
// 	if(mouse_mode_ == TRANS)
// 	{
// 		MeshViewerWidget::mouseMoveEvent(_event);
// 	}
// 	else
// 	{
// 		if( mouse_mode_ != T2_MODE)
// 		{
// 			if( mouse_mode_ == MOVE )
// 			{
// 				move_point_based_lastVertex( _event->x(), _event->y() );
// 				Mesh::Point P(selectedPoint[0],selectedPoint[1],selectedPoint[2]);
// 				mesh.set_point( mesh.vertex_handle(lastestVertex), P );
// 				updateGL();
// 			}
// 		}
// 		else
// 		{
			
// 		}
		
// 	}
// }

// void InteractiveViewerWidget::mouseReleaseEvent(QMouseEvent *_event)
// {
// 	if(mouse_mode_ == TRANS)
// 	{
// 		MeshViewerWidget::mouseMoveEvent(_event);
// 	}
// 	else
// 	{
// 		if(mouse_mode_ != T2_MODE )
// 		{
// 			if( mouse_mode_ == MOVE )
// 			{
// 				move_point_based_lastVertex( _event->x(), _event->y() );
// 				Mesh::Point P(selectedPoint[0],selectedPoint[1],selectedPoint[2]);
// 				mesh.set_point( mesh.vertex_handle(lastestVertex), P );
// 				selectedVertex.clear();
// 				updateGL();
// 			}
// 		}
// 		else
// 		{
// 		}
// 	}
	
// }

// void InteractiveViewerWidget::wheelEvent(QWheelEvent* _event)
// {
// 	if(mouse_mode_ != N_MODE && mouse_mode_ != T2_MODE)
// 	{
// 		MeshViewerWidget::wheelEvent(_event);
// 	}
// }

// void InteractiveViewerWidget::dragEnterEvent(QDragEnterEvent* event)
// {
// 	if( event->mimeData()->hasFormat("text/uri-list") )
// 	{
// 		event->acceptProposedAction();
// 	}
// }

// void InteractiveViewerWidget::dropEvent(QDropEvent* event)
// {
// 	QList<QUrl> urls = event->mimeData()->urls();
// 	if( urls.isEmpty() )
// 		return;
// 	QString fileName = urls.first().toLocalFile();
// 	if (fileName.isEmpty())
// 		return;

// 	if( fileName.endsWith(".off") || fileName.endsWith(".obj") || fileName.endsWith(".stl") || fileName.endsWith(".ply"))
// 	{
// 		if( openMesh(fileName.toLocal8Bit()))
// 		{
// 			emit(loadMeshOK(true,fileName));
// 			setDrawMode(FLAT_POINTS);
// 			setMouseMode(TRANS);

// 			load_parameterization();
// 			input_distortion = para_cutting->get_distortion();
// 			updateGL();
// 		}
// 		else
// 		{
// 			emit(loadMeshOK(false,"No Mesh"));
// 		}
// 	}
// }

 void NInteractiveViewerWidget::load_parameterization(bool silence /*= false*/)
 {
 	if (para_cutting == nullptr)
 	{
 		para_cutting = std::make_unique<ParaQuadCutting>(mesh, str_pathname.c_str(), str_filepath.c_str());
 		M_PE = para_cutting->calc_distortion(silence);

 		selectedFace.clear();
 		std::set<int> selected_v;
 		for (int f : para_cutting->get_flipped_faces())
 		{
 			selectedFace.push_back(f);

 			for (auto fv_h : mesh.fv_range(mesh.face_handle(f))) selected_v.insert(fv_h.idx());
 		}

 		selectedVertex = std::vector<int>(selected_v.begin(), selected_v.end());
 	}
 	else
 	{
 		para_cutting->toggle_viewer_window();
 	}

// 	updateGL();
 }

// void InteractiveViewerWidget::save_parameterization()
// {
// 	if (para_cutting == nullptr) return;

// 	QString fileName = QFileDialog::getSaveFileName(this, tr("Save Parameterization"), QString(para_cutting->get_path().c_str()), tr("OBJ Files (*.obj);;"));
// 	if (!fileName.isEmpty())
// 	{
// 		para_cutting->save_para(fileName.toStdString().c_str());
// 	}

// 	updateGL();
// }

// void InteractiveViewerWidget::reload_mesh()
// {
// 	if (mesh.vertices_empty()) return;

// 	openMesh(str_filepath.toLocal8Bit());
// 	load_parameterization(true);

// 	updateGL();
// }

 void NInteractiveViewerWidget::mirror_neg_charts()
 {
 	if (para_cutting == nullptr) return;

 	para_cutting->flip_neg_charts();
 	para_cutting->get_textured_mesh(mesh);
 	init_properties();

 	para_cutting.reset(nullptr);
 	poly_info.reset(nullptr);
 	initMesh();

 	load_parameterization(true);
	
 	/*updateGL();*/
 }

// void InteractiveViewerWidget::select_long_edges()
// {
// 	if (para_cutting == nullptr) return;

// 	selectedEdge.clear();

// 	for (auto e_h : mesh.edges())
// 	{
// 		if (mesh.is_boundary(e_h)) continue;
		
// 		auto v0 = mesh.to_vertex_handle(mesh.halfedge_handle(e_h, 0));
// 		auto v1 = mesh.to_vertex_handle(mesh.halfedge_handle(e_h, 1));

// //		if (!mesh.is_boundary(v0) && !mesh.is_boundary(v1)) continue;

// 		if (mesh.calc_edge_length(e_h) > 2.0 * mesh_avg_length)
// 		{
// 			selectedEdge.push_back(e_h.idx());
// 		}
// 	}

// 	std::cout << "Select " << selectedEdge.size() << " Edges." << std::endl;

// 	updateGL();
// }

 void NInteractiveViewerWidget::select_cut_bridges()
 {
 	if (para_cutting == nullptr) return;

 	selectedEdge.clear();
 	std::vector<bool> v_oncut(mesh.n_vertices(), false);

 	for (auto e_h : mesh.edges())
 	{
 		if (mesh.property(e_oncut, e_h))
 		{
 			v_oncut[mesh.to_vertex_handle(mesh.halfedge_handle(e_h, 0)).idx()] = true;
 			v_oncut[mesh.to_vertex_handle(mesh.halfedge_handle(e_h, 1)).idx()] = true;
 		}
 	}

 	for (auto e_h : mesh.edges())
 	{
 		if (!mesh.property(e_oncut, e_h))
 		{
 			if (v_oncut[mesh.to_vertex_handle(mesh.halfedge_handle(e_h, 0)).idx()] && v_oncut[mesh.to_vertex_handle(mesh.halfedge_handle(e_h, 1)).idx()])
 			{
 				selectedEdge.push_back(e_h.idx());
 			}
 		}
 	}

 //	std::cout << "Select " << selectedEdge.size() << " Edges." << std::endl;

 	//updateGL();
 }

// void InteractiveViewerWidget::cut_along_seleted()
// {
// 	if (para_cutting == nullptr) return;

// 	if (selectedEdge.empty()) return;

// 	OpenMesh::MPropHandleT<std::vector<Mesh::TexCoord2D>> mvt_list;
// 	OpenMesh::HPropHandleT<int> hvt_index;

// 	if (!mesh.get_property_handle(mvt_list, "mvt_list") || !mesh.get_property_handle(hvt_index, "hvt_index"))
// 	{
// 		std::cout << "Texture data is invalid." << std::endl;
// 		return;
// 	}

// 	std::set<int> new_cut(selectedEdge.begin(), selectedEdge.end());

// 	auto h0 = mesh.halfedge_handle(mesh.edge_handle(selectedEdge[0]), 1);
// 	auto h_iter = h0;
// 	while (!mesh.property(e_oncut, mesh.edge_handle(h_iter)))
// 	{
// 		if (new_cut.count(h_iter.idx() / 2) > 0) h_iter = mesh.prev_halfedge_handle(h_iter);
// 		h_iter = mesh.prev_halfedge_handle(mesh.opposite_halfedge_handle(h_iter));
// 	}

// 	h0 = h_iter;

// 	auto& vt_list = mesh.property(mvt_list);
// 	vt_list.push_back(vt_list[mesh.property(hvt_index, h_iter)]);
// 	mesh.property(hvt_index, h_iter) = vt_list.size() - 1;
// 	while (true)
// 	{
// 		h_iter = mesh.next_halfedge_handle(h_iter);
// 		if (mesh.property(e_oncut, mesh.edge_handle(h_iter)))
// 		{
// 			break;
// 		}
// 		else if (new_cut.count(h_iter.idx() / 2) > 0)
// 		{
// 			vt_list.push_back(vt_list[mesh.property(hvt_index, h_iter)]);
// 			mesh.property(hvt_index, h_iter) = vt_list.size() - 1;
// 			h_iter = mesh.next_halfedge_handle(h_iter);
// 		}
// 		h_iter = mesh.opposite_halfedge_handle(h_iter);
// 		mesh.property(hvt_index, h_iter) = vt_list.size() - 1;
// 	}

// 	selectedEdge.clear();
// 	init_properties();

// 	para_cutting.reset(nullptr);
// 	poly_info.reset(nullptr);
// 	initMesh();

// 	load_parameterization(true);

// 	updateGL();
// }

 void NInteractiveViewerWidget::spilt_seleted_edges()
 {
 	if (para_cutting == nullptr) return;

 	if (selectedEdge.empty()) return;

 	para_cutting->split_edges(selectedEdge);
 	para_cutting->get_textured_mesh(mesh);
 	selectedEdge.clear();
 	init_properties();

 	para_cutting.reset(nullptr);
 	poly_info.reset(nullptr);
 	initMesh();

 	load_parameterization(true);

 	//updateGL();
 }

void NInteractiveViewerWidget::polysquare()
{
	if (para_cutting == nullptr || poly_info != nullptr) return;

	double goal_length = (para_cutting->get_BB_Max() - para_cutting->get_BB_Min()).norm() / 3000.0;
//	double goal_length = mesh_avg_length * 0.03;
	
	std::cout << "--------------------------------------------------\nGoal " << goal_length << std::endl;
	current_phase = 0;

	poly_info = std::make_unique<PolySquareDeformation>(mesh, para_cutting->get_origin_para(), 
		para_cutting->viewer->boundary_vk);
	
	poly_info->set_path(str_pathname.c_str());
	
	if (!poly_info->calc(kernel_width, goal_length, energy_exp_factor, my_file_name))
	{
		poly_info.reset(nullptr);
		return;
	}

	int n_phase = poly_info->get_n_phase();

	//para_viewer_titles.clear();
	//para_viewer_titles.reserve(n_phase);
	//para_viewer_titles.emplace_back("Initialization");
	//para_viewer_titles.emplace_back("Rotation");
	//for (int i = 0; i < (n_phase - 5) / 2; i++)
	//{
	//	para_viewer_titles.emplace_back("Align Deformation " + QString::number(i + 1));
	//	para_viewer_titles.emplace_back("Inner Deformation " + QString::number(i + 1));
	//}
	//para_viewer_titles.emplace_back("Flattening");
	//para_viewer_titles.emplace_back("Untangle");
	//para_viewer_titles.emplace_back("Final");

	change_phase(n_phase - 1);

	/*std::string result_name = str_filename.toStdString().substr(0, str_filename.toStdString().find_last_of('.')) + "_test0.obj";
	std::string result_file = ".\\result2\\" + result_name;
	saveMesh(result_file.c_str());*/
}

void NInteractiveViewerWidget::quad_cutting()
{
	if (para_cutting == nullptr) return;

	// FIXME: load the quad cutter someway
	
	//QString qt_path = QString::fromStdString(para_cutting->get_path());
	//QString qt_path_tri = qt_path + tr("/tri_comp");
	//QString qt_path_quad = qt_path + tr("/quad_comp");

	//QDir dir(qt_path);
// 	if (!QDir(qt_path_tri).exists()) dir.mkdir(qt_path_tri);
// 	if (!QDir(qt_path_quad).exists()) dir.mkdir(qt_path_quad);

	para_cutting->my_cutting(split_thres);
	/*para_cutting->cutting(split_thres);*/

	para_cutting->trans_textured(mesh);

	init_properties();

	para_cutting.reset(nullptr);
	poly_info.reset(nullptr);
	initMesh();

	std::cout << "--------------------------------------------------" << std::endl;
	M_PE = 0.8;
	load_parameterization(false);

	/*std::string result_name = str_filename.toStdString().substr(0, str_filename.toStdString().find_last_of('.')) + "_test2.obj";
	std::string result_file = ".\\result2\\" + result_name;
	saveMesh(result_file.c_str());*/

	//updateGL();
}

void NInteractiveViewerWidget::distortion_reduce()
{
	if (para_cutting == nullptr) return;

	std::cout << "--------------------------------------------------" << std::endl;

	Eigen::MatrixXd v_pos;
	Eigen::MatrixXd uv_v_pos;
	Eigen::MatrixXi fv_id;
	Eigen::MatrixXi uv_fv_id;
	Eigen::MatrixXd f_n;
	para_cutting->get_scaf_info(v_pos, uv_v_pos, fv_id, uv_fv_id, f_n);

	StateManager s_;

	if (M_PE < pe_bound)
	{
		M_PE = std::lround(M_PE * 100) - 1;
		pe_bound = M_PE / 100;
	}
	
	Eigen::MatrixXd v_pos0 = v_pos;
	Eigen::MatrixXd uv_v_pos0 = uv_v_pos;
	Eigen::MatrixXi fv_id0 = fv_id;
	Eigen::MatrixXi uv_fv_id0 = uv_fv_id;
	Eigen::MatrixXd f_n0 = f_n;
	s_.run_interface(v_pos, uv_v_pos, fv_id, uv_fv_id, f_n, input_distortion * 4.0, chart_gap, pe_bound);
	if (s_.is_back)
	{
		StateManager s1_;
		s1_.is_back = true;
		s1_.conv_type = true;
		s1_.run_interface(v_pos0, uv_v_pos0, fv_id0, uv_fv_id0, f_n0, input_distortion * 4.0, chart_gap, pe_bound);
		para_cutting->load_from_scaf(s1_.get_uv());
	}
	else
	{
		para_cutting->load_from_scaf(s_.get_uv());
	}
	//para_cutting->load_from_scaf(s_.get_uv());
	pe_bound = 0.8;
	para_cutting->update_textured_mesh(mesh, true);
	init_properties();

	//std::cout << "str_pathname===" << str_pathname.toStdString() << std::endl;
	//std::cout << "str_filepath===" << str_filepath.toStdString() << std::endl;
	//std::cout << "str_filename===" << str_filename.toStdString() << std::endl;

	//std::string result_name = str_filename.toStdString().substr(0, str_filename.toStdString().find_last_of('.')) + "_result.obj";
	//std::string result_file = ".\\result_bad\\" + result_name;
	//std::string result_file = ".\\result2\\" + result_name;
	//saveMesh(result_file.c_str());

	para_cutting.reset(nullptr);
	poly_info.reset(nullptr);
	initMesh();

	std::cout << "--------------------------------------------------" << std::endl;
	load_parameterization(false);

	//updateGL();
}

void NInteractiveViewerWidget::pre_process()
{
	mirror_neg_charts();
	select_cut_bridges();
	spilt_seleted_edges();
}

//my

void NInteractiveViewerWidget::my_process()          
{
	std::vector<std::string> files;
	__int64 hFile = 0;
	struct __finddata64_t fileinfo;
	hFile = _findfirst64(".\\input4\\*.obj", &fileinfo);
	while (!_findnext64(hFile, &fileinfo))
	{
		files.push_back(fileinfo.name);
	}
	_findclose(hFile);
	//std::cout << "file=========: "<<files.size() << std::endl;
	//std::cout << files[0] << std::endl;

	for (int i = 0; i < files.size(); i++)
	{
		std::string in_file_name = files[i];
		openMesh(in_file_name.c_str());
		//std::cout << files[i] <<std::endl;

		std::string save_path_pre = ".\\result2\\";
		std::string out_file_name = save_path_pre + files[i].substr(0, files[i].find_last_of('.')) + "_axis.obj";
		my_file_name = out_file_name;

		pre_process();
		polysquare();
		quad_cutting();
		distortion_reduce();
	}
}

//

 void NInteractiveViewerWidget::change_phase(int phase)
 {
 	if (poly_info == nullptr) return;

 	int n_phases = para_viewer_titles.size();
 	int new_phase = (phase + n_phases) % n_phases;

 	/*para_cutting->set_viewer_title((para_viewer_titles[new_phase]).c_str());*/
 	poly_info->set_phase(new_phase);

 	int scale_id = n_phases - 4;
 	double factor = 1.0;
 	if (current_phase <= scale_id && new_phase > scale_id) factor = 1.0 / poly_info->get_goal_length();
 	if (current_phase > scale_id && new_phase <= scale_id) factor = poly_info->get_goal_length();
 	para_cutting->update_para(factor);
 	para_cutting->update_textured_mesh(mesh);

 	current_phase = new_phase;
 	//updateGL();
 }

// void InteractiveViewerWidget::dijkstra_path(int src, int dst)
// {
// 	std::set<int> v_cut;

// 	for (auto e_h : mesh.edges())
// 	{
// 		if (!mesh.property(e_oncut, e_h)) continue;

// 		auto v0 = mesh.to_vertex_handle(mesh.halfedge_handle(e_h, 0));
// 		auto v1 = mesh.to_vertex_handle(mesh.halfedge_handle(e_h, 1));

// 		v_cut.insert(v0.idx());
// 		v_cut.insert(v1.idx());
// 	}

// 	Mesh& para = para_cutting->get_origin_para();
// 	OpenMesh::MPropHandleT<std::vector<Mesh::TexCoord2D>> mvt_list;
// 	mesh.get_property_handle(mvt_list, "mvt_list");
// 	auto get_point = [&](int v_id)
// 	{
// 		return para.point(para.vertex_handle(v_id));
// 	};

// 	auto v_src = mesh.vertex_handle(src);

// 	std::map<int, std::pair<double, int>> v_dist;
// 	std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> dist_heap;

// 	dist_heap.emplace(0.0, src);
// 	v_dist[src] = { 0.0, -1 };

// 	while (!dist_heap.empty())
// 	{
// 		auto top = dist_heap.top();
// 		dist_heap.pop();

// 		int v0 = top.second;
// 		double c0 = top.first;

// 		if (c0 > v_dist[v0].first) continue;
		
// 		if (v0 == dst) break;

// 		for (auto vih : mesh.vih_range(mesh.vertex_handle(v0)))
// 		{
// 			int v1 = mesh.from_vertex_handle(vih).idx();

// 			if (mesh.property(e_oncut, mesh.edge_handle(vih)) || (v_cut.count(v1) == 1 && v1 != dst)) continue;

// 			double dc = (get_point(v0) - get_point(v1)).norm();

// 			if (v_dist.find(v1) == v_dist.end() || dc + v_dist[v0].first < v_dist[v1].first)
// 			{
// 				v_dist[v1] = { dc + v_dist[v0].first, vih.idx() };
// 				dist_heap.emplace(v_dist[v1].first, v1);
// 			}
// 		}
// 	}

// 	if (v_dist.find(dst) == v_dist.end())
// 	{
// 		std::cout << src << "-" << dst << " cannot reach. " << std::endl;
// 	}
// 	else
// 	{
// 		OpenMesh::Vec3d v_straight = get_point(dst) - get_point(src);
// 		int v_tag = (std::abs(v_straight[0]) < std::abs(v_straight[1])) ? 0 : 1;
// 		v_straight[v_tag] = 0.0;

// 		v_straight /= v_dist[dst].first;

// 		int v_cur = dst;
// 		while (v_dist[v_cur].second != -1)
// 		{
// 			auto& p = para.point(para.vertex_handle(v_cur));
// 			p[0] = get_point(src)[0] + v_dist[v_cur].first * v_straight[0];
// 			p[1] = get_point(src)[1] + v_dist[v_cur].first * v_straight[1];

// 			selectedEdge.push_back(v_dist[v_cur].second / 2);
// 			mesh.property(e_oncut, mesh.edge_handle(v_dist[v_cur].second / 2)) = true;
// 			v_cur = mesh.to_vertex_handle(mesh.halfedge_handle(v_dist[v_cur].second)).idx();
// 		}
// 	}
// }

// void InteractiveViewerWidget::pick_vertex(int x,int y)
// {
// 	int r = find_vertex_using_selected_point();
// 	lastestVertex = r;
// 	std::cout << "Select Vertex : " << r << " (" << mesh.point(mesh.vertex_handle(r)) << ")" << std::endl;
// 	std::vector<int>::iterator it;
// 	if( (it = std::find(selectedVertex.begin(),selectedVertex.end(), r)) == selectedVertex.end() )
// 	{
// 		selectedVertex.push_back(r);
// 	}
// 	else
// 	{
// 		selectedVertex.erase(it);
// 	}
// 	updateGL();
// }
// void InteractiveViewerWidget::pick_face(int x,int y)
// {
// 	int desiredFace = find_face_using_selected_point();
// 	if(desiredFace < 0) return;
// 	lastestFace = desiredFace;
// 	std::cout << "Select Face : " << desiredFace << " (" << mesh.normal(mesh.face_handle(desiredFace)) << ")" << std::endl;
// 	std::vector<int>::iterator it;
// 	if( (it = std::find(selectedFace.begin(),selectedFace.end(),desiredFace)) == selectedFace.end() )
// 	{
// 		selectedFace.push_back(desiredFace);
// 	}
// 	else
// 	{
// 		selectedFace.erase(it);
// 	}
// 	updateGL();
// }
// void InteractiveViewerWidget::pick_edge(int x,int y)
// {
// 	int desiredEdge = find_edge_using_selected_point();
// 	if(desiredEdge < 0) return;
// 	lastestEdge = desiredEdge;

// 	auto e_h = mesh.edge_handle(desiredEdge);
// 	auto h0 = mesh.halfedge_handle(e_h, 0);
// 	auto h1 = mesh.halfedge_handle(e_h, 1);
// 	QString e_str;
// 	e_str.sprintf("Select Edge : %d (0:%d,%d 1:%d,%d)", desiredEdge, h0.idx(), mesh.face_handle(h0).idx(), h1.idx(), mesh.face_handle(h1).idx());
// 	std::cout << e_str.toStdString() << std::endl;
// 	std::vector<int>::iterator it;
// 	if( (it = std::find(selectedEdge.begin(),selectedEdge.end(),desiredEdge)) == selectedEdge.end() )
// 	{
// 		selectedEdge.push_back(desiredEdge);
// 	}
// 	else
// 	{
// 		selectedEdge.erase(it);
// 	}
// 	updateGL();
// }
// void InteractiveViewerWidget::pick_point(int x,int y)
// {
// 	GLint viewport[4];
// 	glGetIntegerv(GL_VIEWPORT, viewport);
// 	GLdouble winX = double(x);
// 	GLdouble winY = double( height() - y );
// 	GLfloat winZ = 0.0;
// 	glReadPixels((int)winX, (int)winY, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);
// 	gluUnProject(winX, winY, (GLdouble)winZ, &ModelViewMatrix[0], &ProjectionMatrix[0], viewport, &selectedPoint[0], &selectedPoint[1], &selectedPoint[2]);
// }

// void InteractiveViewerWidget::move_point_based_lastVertex(int x,int y)
// {
// 	if(lastestVertex<0 || lastestVertex>=mesh.n_vertices())
// 	{
// 		return;
// 	}
// 	GLint viewport[4];
// 	glGetIntegerv(GL_VIEWPORT, viewport);
// 	GLdouble winX = 0.0;
// 	GLdouble winY = 0.0;
// 	GLdouble winZ = 0.0;
// 	OpenMesh::Vec3d p = mesh.point(mesh.vertex_handle(lastestVertex));
// 	gluProject(p[0], p[1], p[2],  &ModelViewMatrix[0], &ProjectionMatrix[0], viewport, &winX, &winY, &winZ);
	
// 	gluUnProject((GLdouble)(x), (GLdouble)( height() - y ), winZ,  &ModelViewMatrix[0], &ProjectionMatrix[0], viewport, &selectedPoint[0], &selectedPoint[1], &selectedPoint[2]);
// }

// int InteractiveViewerWidget::find_vertex_using_selected_point()
// {
// 	ANNpoint tp = annAllocPt(3); tp[0] = selectedPoint[0]; tp[1] = selectedPoint[1]; tp[2] = selectedPoint[2];
// 	ANNidxArray nnIdx = new ANNidx[1]; ANNdistArray dists = new ANNdist[1];
// 	kdTree->annkSearch(tp, 1, nnIdx, dists);
// 	return nnIdx[0];
// }

// int InteractiveViewerWidget::find_face_using_selected_point()
// {
// 	int rv = find_vertex_using_selected_point();
// 	Mesh::VertexFaceIter vf_it = mesh.vf_iter( mesh.vertex_handle(rv) );
// 	int desiredFace = -1; //double minLen = 10*radius();
// 	std::vector<OpenMesh::Vec3d> tri_p(3); int tri_count = 0;
// 	Mesh::Point resultP(selectedPoint[0],selectedPoint[1],selectedPoint[2]);
// 	for( vf_it; vf_it; ++vf_it )
// 	{
// 		tri_count = 0;
// 		for(Mesh::FaceVertexIter fv_it = mesh.fv_iter(vf_it.handle()); fv_it; ++fv_it)
// 		{
// 			tri_p[tri_count] = mesh.point(fv_it); ++tri_count;
// 		}
// 		if( check_in_triangle_face(tri_p, resultP) )
// 		{
// 			desiredFace = vf_it.handle().idx(); break;
// 		}
// 	}
// 	if(desiredFace < 0)
// 	{
// 		for(Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
// 		{
// 			tri_count = 0;
// 			for(Mesh::FaceVertexIter fv_it = mesh.fv_iter(f_it.handle()); fv_it; ++fv_it)
// 			{
// 				tri_p[tri_count] = mesh.point(fv_it); ++tri_count;
// 			}
// 			if( check_in_triangle_face(tri_p, resultP) )
// 			{
// 				desiredFace = f_it.handle().idx(); break;
// 			}
// 		}
// 	}

// 	return  desiredFace;
// }

// int InteractiveViewerWidget::find_edge_using_selected_point()
// {
// 	int desiredFace = find_face_using_selected_point(); if(desiredFace < 0) return -1;
// 	Mesh::FaceHandle fh = mesh.face_handle(desiredFace);
// 	double min_len= 1e30; int desiredEdge = -1;
// 	Mesh::Point resultP(selectedPoint[0],selectedPoint[1],selectedPoint[2]);
// 	for(Mesh::FaceHalfedgeIter fhe_it = mesh.fh_iter(fh); fhe_it; ++fhe_it)
// 	{
// 		OpenMesh::Vec3d s = mesh.point( mesh.from_vertex_handle(fhe_it) );
// 		OpenMesh::Vec3d e = mesh.point( mesh.to_vertex_handle(fhe_it) );
// 		double dis = OpenMesh::cross(resultP - s, resultP - e).norm() / (s - e).norm();
// 		if(dis < min_len){ min_len = dis; desiredEdge = mesh.edge_handle(fhe_it.handle()).idx(); }
// 	}
	
// 	return desiredEdge;
// }

// void InteractiveViewerWidget::buildIndex()
// {
// 	if(mesh.n_vertices() == 0)
// 		return;

// 	Mesh::VertexIter v_it(mesh.vertices_begin());
// 	Mesh::VertexIter v_end(mesh.vertices_end());
// 	Mesh::Point p;
// 	unsigned nv = mesh.n_vertices();
// 	ANNpointArray dataPts = annAllocPts(nv, 3);
// 	int count = 0;
// 	for(; v_it != v_end; ++v_it)
// 	{
// 		p = mesh.point(v_it);
// 		dataPts[count][0] = p[0]; dataPts[count][1] = p[1]; dataPts[count][2] = p[2];
// 		++count;
// 	}

// 	if(kdTree) delete kdTree;
// 	kdTree = new ANNkd_tree(dataPts, nv, 3);
// }

// //with the first mesh
// void InteractiveViewerWidget::draw_interactive_portion(int drawmode)
// {
// 	glViewport ( 0,0, width(),height());
// 	glMatrixMode( GL_PROJECTION );
// 	glLoadMatrixd( &ProjectionMatrix[0] );
// 	glMatrixMode( GL_MODELVIEW );
// 	glLoadMatrixd( &ModelViewMatrix[0] );

	
// 	{
// 		//draw select vertex, face, edge.
// 		glDisable(GL_LIGHTING);
// 		glDisable(GL_TEXTURE_2D);

// 		glPointSize(1);

// 		switch(mouse_mode_)
// 		{
// 		case POINTPICK:
// 			draw_selected_point();
// 			break;
// 		case VERTEXPICK:
// 			draw_selected_vertex();
// 			break;
// 		case FACEPICK:
// 			draw_selected_face();
// 			break;
// 		case EDGEPICK:
// 			draw_selected_edge();
// 			break;
// 		default:
// 			draw_selected_vertex();
// 			draw_selected_face();
// 			draw_selected_edge();
// 			break;
// 		}
// 	}

// 	if(draw_new_mesh)
// 	{
// 		draw_scene_mesh(drawmode);
// 	}
// }

// //with the second mesh
// void InteractiveViewerWidget::draw_interactive_portion_mesh2()
// {
// 	return;
// }

// void InteractiveViewerWidget::draw_selected_point()
// {
// 	glColor3f(1.0, 0.5, 0.0);
// 	glPointSize(10);
// 	glBegin(GL_POINTS);
// 	glVertex3d(selectedPoint[0],selectedPoint[1],selectedPoint[2]);
// 	glEnd();
// 	glPointSize(1);
// }

// void InteractiveViewerWidget::draw_selected_vertex()
// {
// 	if( selectedVertex.size() > 0 )
// 	{
// 		Mesh::Point p;
// 		glColor3f(1.0, 0.5, 0.0);
// 		glPointSize(12);
// 		glBegin(GL_POINTS);
// 		for(unsigned int i=0;i<selectedVertex.size();++i)
// 		{
// 			p = mesh.point( mesh.vertex_handle(selectedVertex[i]) );
// 			glVertex3dv(p.data());
// 		}
// 		glEnd();
// 		glPointSize(1);
// 	}
// }

// void InteractiveViewerWidget::draw_selected_face()
// {
// 	if( selectedFace.size() > 0 )
// 	{
// 		glColor3f(1.0, 0.5, 1.0);
// 		Mesh::Point p;
// 		Mesh::ConstFaceVertexIter fv_it;
// 		Mesh::FaceHandle f_handle;
// 		for( unsigned int i=0; i<selectedFace.size(); ++i )
// 		{
// 			f_handle = mesh.face_handle(selectedFace[i]);
// 			fv_it = mesh.fv_iter(f_handle);
// 			glBegin(GL_POLYGON);
// 			for( fv_it; fv_it; ++fv_it )
// 			{
// 				glVertex3dv(&mesh.point(fv_it)[0]);
// 			}
// 			glEnd();
// 		}
// 	}
// }

// void InteractiveViewerWidget::draw_selected_edge()
// {
// 	if( selectedEdge.size() > 0)
// 	{
// 		glColor3f(1.0, 0.5, 1.0);
// 		Mesh::Point p1; Mesh::Point p2;
// 		Mesh::EdgeHandle e_handle;
// 		Mesh::HalfedgeHandle he_handle;
// 		for(unsigned int i=0;i<selectedEdge.size();++i)
// 		{
// 			e_handle = mesh.edge_handle(selectedEdge[i]);
// 			he_handle = mesh.halfedge_handle( e_handle, 0 );
// 			p1 = mesh.point( mesh.from_vertex_handle( he_handle ) );
// 			p2 = mesh.point( mesh.to_vertex_handle( he_handle ) );
// 			glBegin(GL_LINES);
// 			glVertex3dv( p1.data() );
// 			glVertex3dv( p2.data() );
// 			glEnd();
// 		}
// 	}
// }

// void InteractiveViewerWidget::draw_scene(int drawmode)
// {
// 	if (!mesh.n_vertices()) { return; }
// 	draw_interactive_portion_mesh2();
// 	draw_interactive_portion(drawmode);

// 	if( !draw_new_mesh )
// 	{
// 		MeshViewerWidget::draw_scene(drawmode);
// 	}
// }

// void InteractiveViewerWidget::render_text_slot(OpenMesh::Vec3d pos, QString str)
// {
// 	/*GLdouble  winX, winY, winZ;
// 	GLint     viewport[4];
// 	glGetIntegerv(GL_VIEWPORT, viewport);
// 	gluProject(pos[0],pos[1],pos[2],&ModelViewMatrix[0][0],&ProjectionMatrix[0][0],viewport,&winX,&winY,&winZ);
// 	int x = (long)winX;
// 	int y = viewport[3]-(long)winY;
// 	render_text(x,y,str);*/
// 	render_text(pos[0],pos[1],pos[2],str);
// }
