#pragma once



#include "common_defines.h"
#include "curvaturecalculator.h"
#include "discretefairer.h"
#include <fstream>
#include <vector>



namespace metrics {



  class CurvatureBased
  {
    struct VertexRecord
    {
      int vh_idx;
      double target_curvature;
      double real_curvature;
    };

    using DataFrame = std::vector<VertexRecord>;
    using SessionData = std::vector<DataFrame>;

    SessionData session_data;
    
    struct DataFrameStatistic
    {
      double absmean_curvature_error;
      double max_curvature_error;
    };
    using SessionStatistic = std::vector<DataFrameStatistic>;

    const core::DiscreteFairer& m_df;

    DataFrameStatistic analizedataFrame(const DataFrame& data_frame) const;
    
  public:
    static size_t global_session_id;

    size_t m_session_id = 0;

    size_t m_iteration_count = 0;

    common::MyMesh& subject; //TODO const?

    std::fstream fw;

    CurvatureBased(common::MyMesh& _subject, const core::DiscreteFairer& df);

    void startSession();

    void postIteration();

    void endSession();



    
  };

  
  
}
