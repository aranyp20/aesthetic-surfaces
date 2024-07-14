#include "curvaturecalculator.h"

#include "common_defines.h"
#include "lsq-plane.hh"
#include <algorithm>
#include <iterator>
#include <regex>


namespace core{
  Eigen::Vector3d tmp_offset;
  size_t dcallnum = 0;

  void CurvatureCalculator::tessellateSurface(const size_t resolution, DerResults Ss) const
  {

	  
#define uvBOUND 2

    const double u_min = -uvBOUND;
    const double u_max = uvBOUND;
    const double v_min = -uvBOUND;
    const double v_max = uvBOUND;

    const double u_step_distance = (u_max - u_min) / resolution;
    const double v_step_distance = (v_max - v_min) / resolution;

    std::vector<Eigen::Vector3d> values;
    for (size_t i = 0; i < resolution; i++)
      {
	const double u = u_min + i * u_step_distance;
	for (size_t j = 0; j < resolution; j++)
	  {
	    const double v = v_min + j * v_step_distance;
	    values.push_back(S(u, v, Ss));
	  }
      }

    // TODO refactor

    std::fstream fw("tessellations/tessellated-" + std::to_string(++dcallnum) + ".obj", std::ios::out);

    fw << "# Vertices\n";

    for (const auto &vertex : values)
      {
	fw << "v " << vertex[0] << " " << vertex[1] << " " << vertex[2] << std::endl;
      }

    fw << "\n# Faces\n";

    for (size_t i = 0; i < (resolution - 1) * (resolution - 1); i++)
      {
	const size_t left_right_index = std::floor(i / (resolution - 1)) + i;

	fw << "f " << left_right_index + 1 << " " << left_right_index + resolution + 1 << " " << left_right_index + resolution + 2 << std::endl;
	fw << "f " << left_right_index + 1 << " " << left_right_index + resolution + 2 << " " << left_right_index + 2 << std::endl;
      }

  }
  
}

namespace core
{
  
    // TODO before integrated: cps.P size check + origo center precondition
    std::vector<CurvatureCalculator::SurfaceParam> CurvatureCalculator::InputPoints::calcUVs(bool use_adaptive_uvs) const
    {
      /*
     std::vector<SurfaceParam> result;

     const auto uvs = LSQPlane::projectToBestFitPlane(P);

     for(const auto& uv : uvs) {
       result.push_back({uv[0], uv[1]});
     }

     return result;
      */
     




      ////////////////////////////////////////////
      
      if (!use_adaptive_uvs) {
	std::vector<SurfaceParam> result;

	for(size_t i=0; i< P.size(); i++) {
	  const auto tt = 2 * M_PI * (static_cast<double>(i)/P.size());
	  result.push_back({std::cos(tt), std::sin(tt)});
	}
	//C-h i
	//C-g

	return result;   

      }
      
      
      
        std::vector<SurfaceParam> retval;

        // retval.push_back({0, 0});

        std::vector<double> alphas;
        double alpha_sum = 0;
        for (size_t i = 0; i < P.size(); i++)
        {
            const auto &p_i = P[i];
            const auto &p_i_next = i == P.size() - 1 ? P[0] : P[i + 1];
            const auto alpha = std::acos(p_i.normalized().dot(p_i_next.normalized()));
            alphas.push_back(alpha);
            alpha_sum += alpha;
        }

        const auto alpha_normalizer = 2 * M_PI / alpha_sum;
        for (auto &alpha : alphas)
        {
            alpha *= alpha_normalizer;
        }
        for (size_t i = 0; i < P.size(); i++)
        {
            const auto h_i = P[i].norm();
            double sum_a = 0;
            
            {
                for (size_t j = 0; j < i; j++)
                {
                    sum_a += alphas[j];
                }
            }
            retval.push_back({h_i * std::cos(sum_a), h_i * std::sin(sum_a)});
        }
        return retval;
    }

    CurvatureCalculator::InputPoints CurvatureCalculator::InputPoints::translatePoints() const
    {
      tmp_offset=center;
        auto clone = *this;

        for (auto &cp : clone.P)
        {
            cp -= clone.center;
        }
        clone.center = {0, 0, 0};

        return clone;
    }

  CurvatureCalculator::DerResults CurvatureCalculator::calcDer(const InputPoints &ipp)
    {
        const auto uvs = ipp.calcUVs(use_adaptive_uvs);


        const auto res_size = 5;
        Eigen::MatrixXd A = Eigen::MatrixXd::Random(uvs.size(), res_size);
        Eigen::MatrixXd c = Eigen::MatrixXd::Random(uvs.size(), 3);

        for (int i = 0; i < uvs.size(); i++)
        {
            const auto &uv = uvs[i];
            A(i, 0) = uv.u;
            A(i, 1) = uv.v;
            A(i, 2) = 0.5 * uv.u * uv.u;
            A(i, 3) = uv.u * uv.v;
            A(i, 4) = 0.5 * uv.v * uv.v;

            c(i, 0) = ipp.P[i][0];
            c(i, 1) = ipp.P[i][1];
            c(i, 2) = ipp.P[i][2];
        }

        Eigen::MatrixXd b = Eigen::MatrixXd::Random(res_size, 3);

        if (uvs.size() >= 5)
        {
            // Eigen uses the least squares method for overdetermined systems by default
            b = A.colPivHouseholderQr().solve(c);

	    if (uvs.size()==6) {
	      //std::cout<<"------------"<<std::endl;
	      //std::cout<<(A.transpose() * A).inverse() * A.transpose()<<std::endl;
	      //std::cout<<"------------"<<std::endl;
	      last_M = (A.transpose() * A).inverse() * A.transpose();

	    }
	}
        else
        {
            // Min norm methods
            b = A.transpose() * (A * A.transpose()).inverse() * c;

	    //last_M = A.transpose() * (A * A.transpose()).inverse();
	}

	//	tessellateSurface(100, {b.row(0), b.row(1), b.row(2), b.row(3), b.row(4)});

        return {b.row(0), b.row(1), b.row(2), b.row(3), b.row(4)};
    }

    Eigen::Vector3d CurvatureCalculator::S(const double u, const double v, const DerResults &Ss) const
    {
        const Eigen::Vector3d S00{0, 0, 0};

        return S00 + u * Ss.Su + v * Ss.Sv + 0.5 * u * u * Ss.Suu + u * v * Ss.Suv + 0.5 * v * v * Ss.Svv + tmp_offset;
    }

    void CurvatureCalculator::calcCurvature(const InputPoints &ipp)
    {
	if (dcallnum == 0) {
	  std::ofstream f("pont.obj");
	  for ( auto& i : ipp.P) {
	    f << "v "<<i[0]<<" "<<i[1]<<" "<<i[2]<<std::endl;
	  }
	   f << "v "<<ipp.center[0]<<" "<<ipp.center[1]<<" "<<ipp.center[2]<<std::endl;
	}
	const auto ip = ipp.translatePoints();


        const auto eq = calcDer(ip);
	//std::cout<<"SuSv: "<<eq.Su<<", "<<eq.Sv<<std::endl;

        fundamental_elements.E = eq.Su.dot(eq.Su);
        fundamental_elements.G = eq.Sv.dot(eq.Sv);
        fundamental_elements.F = eq.Su.dot(eq.Sv);

        normal = eq.Su.cross(eq.Sv).normalized();

        
        fundamental_elements.L = eq.Suu.dot(normal);
        fundamental_elements.M = eq.Suv.dot(normal);
        fundamental_elements.N = eq.Svv.dot(normal);
    }

    double CurvatureCalculator::getMeanCurvature() const
    {
        const auto& fe = fundamental_elements;
        return (fe.L * fe.G - fe.M * fe.F + fe.N * fe.E - fe.M * fe.F) / (2 * (fe.E * fe.G - fe.F * fe.F));
    }

    double CurvatureCalculator::getGaussianCurvature() const
    {
        const auto& fe = fundamental_elements;
        return (fe.N * fe.E - 2 * fe.M * fe.F + fe.L * fe.G) / (fe.E * fe.G - fe.F * fe.F);
    }

    double CurvatureCalculator::getMaxPrincipleCurvature() const
    {
        const auto main_curvature = getMeanCurvature();
        const auto gaussian_curvature = getGaussianCurvature();
        return 0; 
    }

/*
    k2 = B - k1
    k1(B - k1) = A

    k1*B - k1*k1 -A= 0

    a=-1
    b=B
    c=-A
*/

    void CurvatureCalculator::execute(common::MyMesh::VertexHandle &vh)
    {
        const common::MyMesh::Point& vertexPos = mesh.point(vh);

        InputPoints cip;
        cip.center = Eigen::Vector3d(vertexPos[0], vertexPos[1], vertexPos[2]);
	////////////
	for (auto neighborVh : mesh.vv_range(vh)) {

            const auto neighborPos = mesh.point(neighborVh);

            cip.P.emplace_back(neighborPos[0], neighborPos[1], neighborPos[2]);
        }

	//TODO delete
	//std::reverse(cip.P.begin(), cip.P.end());
	

	////////////////////
	/*
	const size_t zone_size = 1;
	std::set<common::MyMesh::VertexHandle> extended_neighborhood;
	std::set<common::MyMesh::VertexHandle> extended_neighborhood2;


	extended_neighborhood.insert(vh);
	
	for (size_t i = 0; i < zone_size; i++) {
	  for (auto e2_neighborVh : extended_neighborhood) {
	    for (auto e_neighborVh : mesh.vv_range(e2_neighborVh)) {
	      extended_neighborhood2.insert(e_neighborVh);
	    }	    
	  }
	  //for(const auto& tt : extended_neighborhood2) {
	  //  extended_neighborhood.insert(tt);
	  //}
	  //extended_neighborhood2.clear();
	  extended_neighborhood = extended_neighborhood2;
	}

	extended_neighborhood.insert(vh);
	
	for (auto collected_neighbor : extended_neighborhood) {
	  const auto collected_neighbor_pos = mesh.point(collected_neighbor);
	  cip.P.emplace_back(collected_neighbor_pos[0], collected_neighbor_pos[1], collected_neighbor_pos[2]);
	}
	*/
	/////////////////////
	
	calcCurvature(cip);
	
        //std::cout<<"curvature: "<<getMeanCurvature()<<" at: "<<dcallnum<<std::endl;
    }

    void CurvatureCalculator::execute(const Eigen::Vector3d vertex_pos, const std::vector<Eigen::Vector3d>& neighbors)
    {
      
        InputPoints cip;
        cip.center = Eigen::Vector3d(vertex_pos[0], vertex_pos[1], vertex_pos[2]);
        for (const auto& neighbor_pos : neighbors)
        {
            cip.P.emplace_back(neighbor_pos[0], neighbor_pos[1], neighbor_pos[2]);
        }

        calcCurvature(cip);
    }


    CurvatureCalculator::FundamentalElements CurvatureCalculator::getFundamentalElements() const
    {
        return fundamental_elements;
    }

    Eigen::Vector3d CurvatureCalculator::getNormal() const
    {
        return normal;
    }

  CurvatureCalculator::CurvatureCalculator(common::MyMesh &mesh, bool _use_adaptive_uvs) : mesh(mesh), use_adaptive_uvs(_use_adaptive_uvs)
    {
    }


  Eigen::Matrix<double, 5, 6> CurvatureCalculator::getLastM() const
  {
    return last_M;
  }
}
