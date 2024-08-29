#include "curvaturebased.h"
#include "common_defines.h"
#include "curvaturecalculator.h"
#include "discretefairer.h"
#include <QtWidgets/qwidget.h>
#include <fstream>
#include <string>



namespace metrics {

  size_t CurvatureBased::global_session_id = 0;

  CurvatureBased::CurvatureBased(common::MyMesh& _subject, const core::DiscreteFairer& df) : subject(_subject), m_df(df)
  {
  }
  

  void CurvatureBased::startSession()
  {
    m_session_id = global_session_id++;
    m_iteration_count = 0;
    session_data.clear();
    
    std::cout << "Metrics session started with id {"<< m_session_id <<"}"<<std::endl;


    

    fw = std::fstream("metrics/curvaturebased-" + std::to_string(m_session_id) + ".txt", std::ios::out);

    fw << "Session {"<<global_session_id<<"}\n"<<std::endl;
    
    
  }

  void CurvatureBased::postIteration()
  {

    DataFrame data_frame;

    core::CurvatureCalculator cc_mirror(subject, true);
    for (common::MyMesh::VertexIter v_it = subject.vertices_begin(); v_it != subject.vertices_end(); ++v_it) {
      common::MyMesh::VertexHandle vh = *v_it;
      if(m_df.extended_vertex_static_infos.at(vh).is_original_vertex) {
	continue;
      }
      VertexRecord vr;
      vr.vh_idx = vh.idx();
      throw std::runtime_error("not implemented");
      //vr.target_curvature = m_df.calcTargetCurvature(m_df.extended_vertex_static_infos.at(vh).weighed_effectors);
      cc_mirror.execute(vh);
      vr.real_curvature = cc_mirror.getMeanCurvature();
      data_frame.push_back(vr);
    }

    session_data.push_back(data_frame);



    m_iteration_count++;
  }

  CurvatureBased::DataFrameStatistic CurvatureBased::analizedataFrame(const DataFrame& data_frame) const
  {
    size_t vert_count = 0;
    double absmean_error = 0;
    double max_abs_error = 0;
    int max_error_vert_idx = 0; 
    
    for(const auto& vert_data : data_frame) {
      const auto vert_abs_error = std::fabs(vert_data.real_curvature - vert_data.target_curvature);
      vert_count++;
      absmean_error += vert_abs_error;
      max_abs_error = std::max(max_abs_error, vert_abs_error);
      max_error_vert_idx = vert_data.vh_idx;
    }
    absmean_error /= vert_count;
    
    DataFrameStatistic retval;
    retval.absmean_curvature_error = absmean_error;
    retval.max_curvature_error = max_abs_error;
    retval.max_curvarure_error_vert_idx = max_error_vert_idx;
    
    return retval;
  }

  void CurvatureBased::endSession()
  {
    for(size_t iteration = 0; iteration < session_data.size(); iteration++) {
      fw << "Iteration {" << iteration << "}\n";
      const auto dfs = analizedataFrame(session_data[iteration]);
      fw << "Absmean error: " << dfs.absmean_curvature_error << "\n";
      fw << "Max error: " << dfs.max_curvature_error << "\n";
      fw << "Worst vertex: " << dfs.max_curvarure_error_vert_idx << "\n";
      fw << "----------------" << std::endl;
    }
  }

  
}
