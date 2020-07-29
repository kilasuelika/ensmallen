/** @author Zhou Yao


*/
#ifndef _EIGEN_HELPER_RAND_
#define _EIGEN_HELPER_RAND_
#include<random>

namespace eigen_helper{
	
	class RandGen{
		std::random_device rd;
	
	public:
		//Random normal distribution.
		template<typename Mat, typename S=double>
		void randn(Mat& m,S mu=0, S var=1){
			typedef typename Mat::Scalar ElemType;
			std::normal_distribution<ElemType> dist{mu,var};
			
			for(int i=0;i<m.rows();++i){
				for(int j=0;j<m.cols();++j){
					m.coeffRef(i,j)=dist(rd);
				};
			};
		};
	};
	
};
#endif