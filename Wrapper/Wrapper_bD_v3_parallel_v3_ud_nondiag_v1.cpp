// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;



// [[Rcpp::export]]


arma::mat upper2vec( arma::mat adj, int include_diag  ){
	
	double p = adj.n_rows;
	double nedge = p*(p-1)/2;
	arma::mat gij;
	
	int index = 0;
	if( include_diag == 0 ){
		gij = arma::zeros(nedge,1);
		for( int j = 1; j < p; j++ ){
			for( int i = 0; i < j; i++ ){
				gij(index) = adj(i,j);
				index++;
			}
		}
	}else{
		gij = arma::zeros(nedge+p,1);
		for( int j = 0; j < p; j++ ){
			for( int i = 0; i <= j; i++ ){
				gij(index) = adj(i,j);
				index++;
			}
		}
	}
	
	return(gij);
	
}




// [[Rcpp::export]]


arma::mat vec2upper( arma::mat vec, int include_diag ){
	
	double r = vec.n_rows;
	double c = vec.n_cols;
	double nedge = std::max(r,c);
	double p;
	
	if( include_diag == 0 ){
		p = 0.5 + 0.5*std::sqrt(1+8*nedge);
	}else{
		p = -0.5 + 0.5*std::sqrt(1+8*nedge);
		nedge = nedge - p;
	}
	
	arma::mat Beta = arma::zeros(p,p);
	int index = 0;
	
	if( include_diag == 0 ){
		for( int j = 1; j < p; j++ ){
			for( int i = 0; i < j; i++ ){
				Beta(i,j) = vec(index);
				index++;
			}
		}
	}else{
		for( int j = 0; j < p; j++ ){
			for( int i = 0; i <= j; i++ ){
				Beta(i,j) = vec(index);
				index++;
			}
		}
	}
	
	return(Beta);
	
}


// [[Rcpp::export]]

double log_multi_gamma(double p, double n) {
	
	double f = (p*(p-1)/4)*log(arma::datum::pi);
	
	for( double j = 1; j <= p; j++ ){
		f += lgamma(n+(1-j)/2);
	}
	
	return f;
	
}




// [[Rcpp::export]]

double log_iwishart_InvA_const_parallel(double df, arma::mat Sr) {
	
	double p = Sr.n_rows;
	arma::mat S (Sr.begin(),p,p);
	
	double log_detS = log(arma::det(S));
	
	//double iwc = 0;
	//if( std::isnan(log_detS)==0 ){
	double iwc = (df+p-1)*(log_detS - p*log(2))/2-log_multi_gamma(p,(df+p-1)/2);
	//}
	
	return iwc;
	
}



// [[Rcpp::export]]

double log_J_parallel(double h, arma::mat Br, double a11) {
	
	double p = Br.n_rows;
	//NumericMatrix cppB = clone(Br);
	arma::mat B (Br.begin(),p,p);
	
	double B11 = B(0,0), B12 = B(0,1), B22 = B(1,1);
	arma::mat mB22;
	mB22 = B(1,1);
	
	double y = log(2*arma::datum::pi/B22)/2;
	
	return y - log_iwishart_InvA_const_parallel(h,mB22) + (h-1)*log(a11)/2 - (B11-B12*B12/B22)*a11/2;
	
}




// [[Rcpp::export]]

double log_H_parallel(double b_prior, arma::mat D_priorr, double n, arma::mat Sr, arma::mat Cr, double i, double j) {
	
	double p = Cr.n_rows;
//	arma::mat cppC = Cr;
//	arma::mat cppC1 = Cr;
	arma::mat C (Cr.begin(),p,p);
	arma::mat D_prior (D_priorr.begin(),p,p);
	arma::mat S (Sr.begin(),p,p);
	
	//// (i,j) = 0
	
	arma::mat C0 (Cr.begin(),p,p);
	C0(i,j) = 0;
	C0(j,i) = 0;
	
	arma::mat C_120 = C0.row(j);
	C_120.shed_col(j);
	
	arma::mat C_220 = C0;
	C_220.shed_row(j);
	C_220.shed_col(j);
	
	arma::mat invC_220 = arma::inv(C_220);
	arma::mat c = C_120*invC_220*arma::trans(C_120);
	
	double temp[] = {C(i,i),0,0,c(0,0)};
	arma::mat C0_ij (temp,2,2);
	
	//// (i,j) = 1
	
	arma::mat C_12 = arma::join_cols(C.row(i),C.row(j));
	C_12.shed_col(i);
	C_12.shed_col(j-1);
	
	arma::mat C_22 (Cr.begin(),p,p);
	C_22.shed_row(i);
	C_22.shed_row(j-1);
	C_22.shed_col(i);
	C_22.shed_col(j-1);
	arma::mat invC_22 = arma::inv(C_22);
	arma::mat Ce = C_12*invC_22*arma::trans(C_12);
	
	arma::mat C_row = arma::join_cols(C.row(i),C.row(j));
	arma::mat Cee = arma::join_rows(C_row.col(i),C_row.col(j));
	arma::mat A = Cee - Ce;
	
	double a11 = A(0,0);
	
	double b_post = b_prior + n;
	arma::mat D_post = D_prior + S;
	arma::mat D_post_row = arma::join_cols(D_post.row(i),D_post.row(j));
	arma::mat D_postee = arma::join_rows(D_post_row.col(i),D_post_row.col(j));
	arma::mat D_postjj;
	D_postjj = D_post(j,j);
	
	double h1 = -log_iwishart_InvA_const_parallel(b_post,D_postjj);
	double h2 = -log_J_parallel(b_post,D_postee,a11);
	double h3 = (n+b_prior-2)/2*(log(a11));
	double h4 = -trace(D_postee*(C0_ij-Ce))/2;
	
	return(h1+h2+h3+h4);
	
}




// [[Rcpp::export]]

std::map< std::string, std::map<int, std::map< std::string, std::vector<int> > > > makedecompgraph_parallel(arma::mat A){
	
	double p = A.n_rows;
	
//	NumericMatrix A_fix = clone(A);
	arma::mat Adj (A.begin(),p,p);
	
//	NumericMatrix A_fix2 = clone(A);
	arma::mat Adj_temp (A.begin(),p,p);
	arma::vec order = arma::linspace<arma::vec>(0, p-1, p);
	
	//arma::mat a0colsum = arma::zeros(1,p);
	
	//// Corresponding MATLAB code
	// i=1;
	// while i<p
	// 	[a,b] = max(sum(Adj(1:i,i+1:end),1)); 
	// 	order([i+1 b+i]) = order([b+i i+1]);
	// 	i=i+1;
	// 	Adj = full(A); Adj=Adj(order,order);
	// end
	
	int i = 0;
	while( i < p-1 ){
		
		arma::mat a0colsum;
		
		if( i == 0 ){ 
			a0colsum = Adj_temp(0,arma::span(1,p-1)); 
		}else{ 
			a0colsum = arma::sum(Adj_temp(arma::span(0,i),arma::span(i+1,p-1)),0);
		}
		
		int b = std::distance(a0colsum.begin(), std::max_element(a0colsum.begin(),a0colsum.end()));
		
		int orderi1 = order[i+1];
		int orderbi = order[b+i+1];
		order[i+1] = orderbi;
		order[b+i+1] = orderi1;
		
		arma::mat Adj_new0 = Adj.row(order[0]);
		for( int r = 1; r < p; r++ ){
			Adj_new0 = arma::join_cols(Adj_new0,Adj.row(order[r]));
		}
		
		arma::mat Adj_new1 = Adj_new0.col(order[0]);
		for( int c = 1; c < p; c++ ){
			Adj_new1 = arma::join_rows(Adj_new1,Adj_new0.col(order[c]));
		}
		
		Adj_temp = Adj_new1;
		
		i++;
		
	}
	
	

	//// Corresponding MATLAB code
	//
	// numberofcliques = 1; C(1).ID = [1]; i=2;
	// while i<=p
	// 	if(sum(Adj(i,C(numberofcliques).ID))==length(C(numberofcliques).ID))
	// 	  C(numberofcliques).ID = [C(numberofcliques).ID i];
	// 	else
	// 		C(numberofcliques).dim = length(C(numberofcliques).ID);
	// 		numberofcliques = numberofcliques + 1;
	// 		C(numberofcliques).ID = union(i,find(Adj(i,1:i)==1));
	// 	end
	// 	i=i+1;
	// end
	// C(numberofcliques).dim = length(C(numberofcliques).ID);
	
	
	// C(1).ID = [1];
	
	std::map<int, std::map< std::string, std::vector<int> > > C;
	std::vector<int> IDvect(1,0);
	C[0]["ID"] = IDvect;
	
	
	// numberofcliques = 1;
	// i = 2;
	
	int numberofcliques = 0;	// (# of cliques)-1 for indexing
	i = 1;
	
	
	// while i<=p
	// 	if(sum(Adj(i,C(numberofcliques).ID))==length(C(numberofcliques).ID))
	// 	  C(numberofcliques).ID = [C(numberofcliques).ID i];
	// 	else
	// 		C(numberofcliques).dim = length(C(numberofcliques).ID);
	// 		numberofcliques = numberofcliques + 1;
	// 		C(numberofcliques).ID = union(i,find(Adj(i,1:i)==1));
	// 	end
	// 	i=i+1;
	// end
	// C(numberofcliques).dim = length(C(numberofcliques).ID);
	
	while( i < p ){
		
		IDvect = C[numberofcliques]["ID"];
		arma::mat Adjrow = Adj_temp.row(i);
		arma::mat AdjC = Adjrow.col(IDvect[0]);
		std::vector<int> dimvect(1,IDvect.size());
		if( IDvect.size() > 1 ){
			for( int a = 1; a < dimvect[0]; a++ ){
				AdjC = arma::join_cols(AdjC,Adjrow.col(IDvect[a]));
			}
		}
		
		arma::mat sumAdjC = arma::sum(AdjC);
		
		if( sumAdjC(0,0)==IDvect.size() ){
			C[numberofcliques]["ID"].push_back(i);
		}else{
			
			C[numberofcliques]["dim"] = dimvect;
			numberofcliques++;
			
			// union(i,find(Adj(i,1:i)==1));
			arma::uvec oneindex = find(Adj(i,arma::span(0,i))==1);
			arma::uvec nonduplicate = find(oneindex != i);
			
			std::vector<int> tempIDvect(1,i);
			if( nonduplicate.n_elem > 0 ){
				for( int s = 0; s < nonduplicate.n_elem; s++ ){
					tempIDvect.push_back(oneindex(nonduplicate(s)));
				}
			}
			
			if( tempIDvect.size() > 1 ){
				std::sort(tempIDvect.begin(),tempIDvect.end());
			}
			
			C[numberofcliques]["ID"] = tempIDvect;
			
		}
		
		i++;
		
		if( i == p ){
			C[numberofcliques]["dim"] = dimvect;
		}
		
	}
	
	
	//// Corresponding MATLAB code
	//
	// for i=1:numberofcliques
	// 	C(i).ID = sort(order(C(i).ID));
	// 	C(i).names = nodenames(C(i).ID,:);
	// end
	
	std::map<int, std::map< std::string, std::vector<int> > > C_sort;
	
	for( int list = 0; list <= numberofcliques; list++ ){
		C_sort[list] = C[order[list]];
	}
	
	//// Corresponding MATLAB code
	//
	// UN = C(1).ID; S(1).ID=[]; S(1).names=[]; S(1).dim=[];
	//
	// for i=2:numberofcliques
	// 	S(i).ID    = intersect(UN,C(i).ID);
	// 	S(i).dim   = length(S(i).ID);
	// 	S(i).names = nodenames(S(i).ID,:);
	// 	S(i).names = [];
	// 	UN = union(UN,C(i).ID);
	// end
	// C(1).names = nodenames(C(1).ID,:);  S(1).names = []; 
	// G{1}=C; G{2}=S;
	
	std::map<int, std::map< std::string, std::vector<int> > > S;
	
	std::vector<int> UN(C[0]["ID"]);
	S[0]["ID"] = std::vector<int>();
	S[0]["dim"] = std::vector<int>();
	
	for( int c = 1; c <= numberofcliques; c++){
			
		std::vector<int> temp_CID = C_sort[c]["ID"];
		std::vector<int> inter;
		std::set_intersection(UN.begin(),UN.end(),temp_CID.begin(),temp_CID.end(),back_inserter(inter));
		
		S[c]["ID"] = inter;
		S[c]["dim"] = std::vector<int>(1,inter.size());
		
		std::vector<int> uni;
		std::set_union(UN.begin(),UN.end(),temp_CID.begin(),temp_CID.end(),back_inserter(uni));
		UN.erase(UN.begin(),UN.end());
		UN.insert(UN.begin(),uni.begin(),uni.end());
		
	}
	
//	List G;
//	G["C"] = wrap(C_sort);
//	G["S"] = wrap(S);
	
	std::map< std::string, std::map<int, std::map< std::string, std::vector<int> > > > G;
	G["C"] = C_sort;
	G["S"] = S;
	
	return G;
	
}




// [[Rcpp::export]]

double log_hiwishart_InvA_const_parallel( arma::mat adjr, double df, arma::mat Sr) {
	
	double p = adjr.n_rows;
	arma::mat adj (adjr.begin(),p,p);
	arma::mat S (Sr.begin(),p,p);
	
	std::map< std::string, std::map<int, std::map< std::string, std::vector<int> > > > G = makedecompgraph_parallel(adj);
	
	std::map<int, std::map< std::string, std::vector<int> > > cliques = G["C"];
	std::map<int, std::map< std::string, std::vector<int> > > separators = G["S"];
	int numberofcliques = cliques.size();
	
	// S(cliques(1).ID,cliques(1).ID)
	std::map< std::string, std::vector<int> > temp_list = cliques[0];
	std::vector<int> temp_ID = temp_list["ID"];
	
	arma::mat submat_row = S.row(temp_ID[0]);
	for( int i = 1; i < temp_ID.size(); i++ ){
		submat_row = arma::join_cols(submat_row,S.row(temp_ID[i]));
	}
	arma::mat submat = submat_row.col(temp_ID[0]);
	for( int j = 1; j < temp_ID.size(); j++ ){
		submat = arma::join_rows(submat,submat_row.col(temp_ID[j]));
	}
	
	double f = log_iwishart_InvA_const_parallel(df,submat);
	
	//Rcout << "numberofcliques is " << numberofcliques << std::endl;
	
	if( numberofcliques > 1 ){
		for( int c = 1; c < numberofcliques; c++ ){
			
			// S(cliques(i).ID,cliques(i).ID)
			std::map< std::string, std::vector<int> > temp_listC = cliques[c];
			std::vector<int> temp_IDC = temp_listC["ID"];
			
			arma::mat submat_rowC;
			arma::mat submatC;
			
			if( temp_IDC.size() > 0 ){
				submat_rowC = S.row(temp_IDC[0]);
				for( int i = 1; i < temp_IDC.size(); i++ ){
					submat_rowC = arma::join_cols(submat_rowC,S.row(temp_IDC[i]));
				}
				submatC = submat_rowC.col(temp_IDC[0]);
				for( int j = 1; j < temp_IDC.size(); j++ ){
					submatC = arma::join_rows(submatC,submat_rowC.col(temp_IDC[j]));
				}
			}
			
			// S(separators(i).ID,separators(i).ID)
			std::map< std::string, std::vector<int> > temp_listS = separators[c];
			std::vector<int> temp_IDS = temp_listS["ID"];
			
			arma::mat submat_rowS;
			arma::mat submatS;
			
			if( temp_IDS.size() > 0 ){
				submat_rowS = S.row(temp_IDS[0]);
				for( int i = 1; i < temp_IDS.size(); i++ ){
					submat_rowS = arma::join_cols(submat_rowS,S.row(temp_IDS[i]));
				}
				submatS = submat_rowS.col(temp_IDS[0]);
				for( int j = 1; j < temp_IDS.size(); j++ ){
					submatS = arma::join_rows(submatS,submat_rowS.col(temp_IDS[j]));
				}
			}
			
			double f0 = log_iwishart_InvA_const_parallel(df,submatC);
			
			f = f + f0 - log_iwishart_InvA_const_parallel(df,submatS);
			
		}
	}
	
	
	return(f);
	
}




// [[Rcpp::export]]

double log_GWishart_ud_const_mc_serial( double b, arma::mat D, arma::mat adj, int N ){
	
	double p = D.n_rows;
	
	arma::mat T = arma::chol(arma::inv(D));
	arma::mat T1 = T*arma::diagmat(1/T.diag());
	
	arma::mat A = vec2upper(upper2vec(adj,0),0);
	arma::mat nu = arma::trans(sum(A,1));
	
	arma::mat logJeach (1,N,arma::fill::zeros);
	
//	Rcout << "ifelse condition = " << (arma::accu(abs(T1)) - arma::accu(T1.diag())) << std::endl;
	
	if( (arma::accu(abs(T1)) - arma::trace(T1)) != 0.0 ){
		
		for( int iter = 0; iter < N; iter++ ){
			
			arma::mat Psi (p,p,arma::fill::zeros);
			Psi.diag() = arma::sqrt(arma::chi2rnd(b+nu));
			
//			Rcout << "A = \n" << A << std::endl;
//			Rcout << "nu = \n" << nu << std::endl;
//			Rcout << "Psi before update = \n" << Psi << std::endl;
			
			Psi.elem( find(A == 1) ) = arma::randn( arma::accu(nu) ) ;
			
//			Rcout << "Psi after update = \n" << Psi << std::endl;
			
			logJeach(iter) = 0 ;
			
			/////// if i = 1 ///////
			for( double j = 1; j < p; j++ ){
				if( A(0,j) == 0 ){
					Psi(0,j) = -arma::accu( arma::trans(Psi.row(0).cols(0,j-1))%(T1.rows(0,j-1).col(j)) );
					logJeach(iter) = logJeach(iter) - (Psi(0,j)*Psi(0,j))/2;
				}
			}
			
			/////// if i > 1 ///////
			for( double i = 1; i < (p-1); i++ ){
				for( double j = (i+1); j < p; j++ ){
					if( A(i,j) == 0 ){
						Psi(i,j) = -arma::accu( arma::trans(Psi.row(i).cols(i,j-1))%(T1.rows(i,j-1).col(j)) );
						for( double r = 0; r < (i-1); r++ ){
							Psi(i,j) = Psi(i,j) - (1/Psi(i,i))*(Psi(r,i) + arma::accu(arma::trans(Psi.row(r).cols(r,i-1))%T1.rows(r,i-1).col(i)))*(Psi(r,j)+arma::accu(arma::trans(Psi.row(r).cols(r,j-1))%(T1.rows(r,j-1).col(j))));
						}
						logJeach(iter) = logJeach(iter) - (Psi(i,j)*Psi(i,j))/2;
					}
				}
			}
			
		}
		
	}else{
		
		for( int iter = 0; iter < N; iter++ ){
			
			arma::mat Psi (p,p,arma::fill::zeros);
			Psi.diag() = arma::chi2rnd(b+nu);
			
//			Rcout << "A = \n" << A << std::endl;
//			Rcout << "nu = \n" << nu << std::endl;
//			Rcout << "Psi before update = \n" << Psi << std::endl;
			
			Psi.elem( find(A == 1) ) = arma::randn( arma::accu(nu) ) ;
			
//			Rcout << "Psi after update = \n" << Psi << std::endl;
			
			logJeach(iter) = 0 ;
			
			/////// if i > 1 ///////
			for( double i = 1; i < (p-1); i++ ){
				for( double j = (i+1); j < p; j++ ){
					if( A(i,j) == 0.0 ){
						Psi(i,j) = -1/Psi(i,i)*arma::accu(Psi.rows(0,i-1).col(i)%Psi.rows(0,i-1).col(j));
						logJeach(iter) = logJeach(iter) - (Psi(i,j)*Psi(i,j))/2;
					}
				}
			}
			
		}
		
	}
	
	
	arma::mat bi = arma::sum(adj,0);
	double logC = arma::accu( nu/2*log(2*arma::datum::pi) + (b+nu)/2*log(2) + lgamma((b+nu)/2) + (b+bi-1)%arma::trans(log(T.diag())) );
	
	double offset = logJeach.max();
	double logJmc = log( arma::mean(arma::mean(arma::exp(logJeach-offset))) ) + offset;
	double log_c = logJmc + logC;
	
	return(log_c);
	
}





// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace RcppParallel;

struct log_GWishart_ud_const_mc_worker : public Worker
{
	// source matrix
	double br;
	const RMatrix<double> Dr;
	const RMatrix<double> adjr;		// (nedge+p) X numfiles
	int Nr;
	
	// destination matrix
	RVector<double> log_c;
	
	// initialize with source and destination
	log_GWishart_ud_const_mc_worker( double br, const NumericMatrix& Dr, 
	const NumericMatrix& adjr, int Nr, NumericVector log_c)
	: br(br), adjr(adjr), Dr(Dr), Nr(Nr), log_c(log_c) {}
	
	
	void operator()(std::size_t begin, std::size_t end) {
		
		double b = br;
		int N = Nr;
		double p = Dr.nrow();
		double numfiles = adjr.ncol();
		double nedge = p*(p-1)/2;
		
		RMatrix<double> cppD = Dr;
		arma::mat D (cppD.begin(),p,p);
		
		RMatrix<double> cppadj = adjr;
		arma::mat adj_all (cppadj.begin(),nedge+p,numfiles);
		
		
		for( std::size_t l = begin; l < end; l++ ){
			
			arma::mat adj_temp = vec2upper(adj_all.col(l),1);
			arma::mat adj = adj_temp + arma::trans(adj_temp) - arma::eye<arma::mat>(p,p);
			
			//////////////////////////////////////////////////////////////////
			arma::mat T = arma::chol(arma::inv(D));
			arma::mat T1 = T*arma::diagmat(1/T.diag());
			
			arma::mat A = vec2upper(upper2vec(adj,0),0);
			arma::mat nu = arma::trans(sum(A,1));
			
			arma::mat logJeach (1,N,arma::fill::zeros);
			
			if( (arma::accu(abs(T1)) - arma::trace(T1)) != 0.0 ){
				
				for( int iter = 0; iter < N; iter++ ){
					
					arma::mat Psi (p,p,arma::fill::zeros);
					Psi.diag() = arma::sqrt(arma::chi2rnd(b+nu));
					
					Psi.elem( find(A == 1) ) = arma::randn( arma::accu(nu) ) ;
					
					logJeach(iter) = 0 ;
					
					/////// if i = 1 ///////
					for( double j = 1; j < p; j++ ){
						if( A(0,j) == 0 ){
							Psi(0,j) = -arma::accu( arma::trans(Psi.row(0).cols(0,j-1))%(T1.rows(0,j-1).col(j)) );
							logJeach(iter) = logJeach(iter) - (Psi(0,j)*Psi(0,j))/2;
						}
					}
					
					/////// if i > 1 ///////
					for( double i = 1; i < (p-1); i++ ){
						for( double j = (i+1); j < p; j++ ){
							if( A(i,j) == 0 ){
								Psi(i,j) = -arma::accu( arma::trans(Psi.row(i).cols(i,j-1))%(T1.rows(i,j-1).col(j)) );
								for( double r = 0; r < (i-1); r++ ){
									Psi(i,j) = Psi(i,j) - (1/Psi(i,i))*(Psi(r,i) + arma::accu(arma::trans(Psi.row(r).cols(r,i-1))%T1.rows(r,i-1).col(i)))*(Psi(r,j)+arma::accu(arma::trans(Psi.row(r).cols(r,j-1))%(T1.rows(r,j-1).col(j))));
								}
								logJeach(iter) = logJeach(iter) - (Psi(i,j)*Psi(i,j))/2;
							}
						}
					}
					
				}
				
			}else{
				
				for( int iter = 0; iter < N; iter++ ){
					
					arma::mat Psi (p,p,arma::fill::zeros);
					Psi.diag() = arma::chi2rnd(b+nu);
					
					Psi.elem( find(A == 1) ) = arma::randn( arma::accu(nu) ) ;
					
					logJeach(iter) = 0 ;
					
					/////// if i > 1 ///////
					for( double i = 1; i < (p-1); i++ ){
						for( double j = (i+1); j < p; j++ ){
							if( A(i,j) == 0.0 ){
								Psi(i,j) = -1/Psi(i,i)*arma::accu(Psi.rows(0,i-1).col(i)%Psi.rows(0,i-1).col(j));
								logJeach(iter) = logJeach(iter) - (Psi(i,j)*Psi(i,j))/2;
							}
						}
					}
					
				}
				
			}
			
			
			arma::mat bi = arma::sum(adj,0);
			double logC = arma::accu( nu/2*log(2*arma::datum::pi) + (b+nu)/2*log(2) + lgamma((b+nu)/2) + (b+bi-1)%arma::trans(log(T.diag())) );
			
			double offset = logJeach.max();
			double logJmc = log( arma::mean(arma::mean(arma::exp(logJeach-offset))) ) + offset;
			log_c[l] = logJmc + logC;
			/////////////////////////////////////////////////////////////////////
			
		}
		
	}
	
};




// [[Rcpp::export]]

NumericVector log_GWishart_ud_const_mc_parallel( double b, NumericMatrix D, NumericMatrix adj, int N ){
	
	// allocate the output matrix
	NumericVector output(adj.ncol());
	
	// (pass input and output matrixes)
	log_GWishart_ud_const_mc_worker log_GWishart_ud_const_mc(b, D, adj, N, output);
	
	
	// call parallelFor to do the work
	parallelFor(0, adj.ncol(), log_GWishart_ud_const_mc);
	
	// return the output matrix
	return output;
	
}







// [[Rcpp::export]]

double log_GWishart_pdf_parallel( arma::mat Kr, double b, arma::mat Dr, arma::mat adjr) {
	
	double p = adjr.n_rows;
	arma::mat adj (adjr.begin(),p,p);
	arma::mat K (Kr.begin(),Kr.n_rows,Kr.n_cols);
	arma::mat D (Dr.begin(),Dr.n_rows,Dr.n_cols);
	
	double log_pdf_unnormalized = (b-2)*log(arma::det(K))/2-trace(D*K)/2;
	
	double log_pdf = log_pdf_unnormalized + log_hiwishart_InvA_const_parallel(adj,b,D);  
	
	return(log_pdf);
	
}



// [[Rcpp::export]]

double log_GWishart_NOij_pdf_parallel( double b_prior, arma::mat D_priorr, arma::mat Cr, double i, double j, double edgeij ) {
	
	double p = D_priorr.n_rows;
//	NumericMatrix cppC = clone(Cr);
//	arma::mat C (cppC.begin(),p,p,false);
	arma::mat C (Cr.begin(),p,p);
	arma::mat D_prior (D_priorr.begin(),p,p);
	
	double l;
	
	if( edgeij == 0 ){
		
		C(i,j) = 0;
		C(j,i) = 0;
		
		arma::mat C_12 = C.row(j);
		C_12.shed_col(j);
		
		arma::mat C_22 = C;
		C_22.shed_row(j);
		C_22.shed_col(j);
		arma::mat invC_22 = arma::inv(C_22);
		
		arma::mat c = C_12*invC_22*arma::trans(C_12);
		
		arma::mat C_new = C;
		C_new(j,j) = c(0,0);
		
		arma::mat Djj;
		Djj = D_prior(j,j);
		
		double log_detC_22 = log(arma::det(C_22));
		double traceDC = trace(D_prior*C_new);
		
		l = -log_iwishart_InvA_const_parallel(b_prior,Djj) + (b_prior-2)*log_detC_22/2 - traceDC/2;
		
	}else{
		
		arma::mat C_12 = arma::join_cols(C.row(i),C.row(j));
		C_12.shed_col(i);
		C_12.shed_col(j-1);
		
		arma::mat C_22 = C;
		C_22.shed_row(i);
		C_22.shed_row(j-1);
		C_22.shed_col(i);
		C_22.shed_col(j-1);
		arma::mat invC_22 = arma::inv(C_22);
		
		arma::mat C_row = arma::join_cols(C.row(i),C.row(j));
		arma::mat C_cliqueid = arma::join_rows(C_row.col(i),C_row.col(j));
		
		arma::mat A = C_cliqueid - C_12*invC_22*arma::trans(C_12);
		A = (A + arma::trans(A))/2;
		
		double log_detC = log(arma::det(C));
		double log_Joint = (b_prior-2)*log_detC/2 - trace(D_prior*C)/2;
		
		arma::mat D_row = arma::join_cols(D_prior.row(i),D_prior.row(j));
		arma::mat D_cliqueid = arma::join_rows(D_row.col(i),D_row.col(j));
		arma::mat ones2 = arma::ones(2,2);
		double logK2by2 = log_GWishart_pdf_parallel(A,b_prior,D_cliqueid,(ones2));
		
		arma::mat V = arma::inv(D_cliqueid);
		arma::mat D_priorii;
		D_priorii = 1/V(0,0);
		arma::mat A11;
		A11 = A(0,0);
		arma::mat one = arma::ones(1,1);
		double logKii = log_GWishart_pdf_parallel(A11,b_prior+1,D_priorii,(one));
		
		l = log_Joint + logKii - logK2by2;
		
	}
	
	return(l);
	
}



//' @title Generate Random Wishart Distribution
//' @description Creates a random wishart distribution when given degrees of freedom and a sigma matrix. 
//' @param df An \code{int}, which gives the degrees of freedom of the Wishart.  (> 0)
//' @param S A \code{matrix} with dimensions m x m that provides Sigma, the covariance matrix. 
//' @return A \code{matrix} that is a Wishart distribution, aka the sample covariance matrix of a Multivariate Normal Distribution
//' @seealso \code{\link{riwishart}} 
//' @author James J Balamuta
//' @examples 
//' #Call with the following data:
//' rwishart(3, diag(2))
//' 
//' # Validation
//' set.seed(1337)
//' S = toeplitz((10:1)/10)
//' n = 10000
//' o = array(dim = c(10,10,n))
//' for(i in 1:n){
//' o[,,i] = rwishart(20, S)
//' }
//' mR = apply(o, 1:2, mean)
//' Va = 20*(S^2 + tcrossprod(diag(S)))
//' vR = apply(o, 1:2, var)
//' stopifnot(all.equal(vR, Va, tolerance = 1/16))
//' 
// [[Rcpp::export]]
arma::mat rwishart(unsigned int df, const arma::mat& S){
  // Dimension of returned wishart
  unsigned int m = S.n_rows;
  
  // Z composition:
  // sqrt chisqs on diagonal
  // random normals below diagonal
  // misc above diagonal
  arma::mat Z(m,m);
  
  // Fill the diagonal
  for(unsigned int i = 0; i < m; i++){
    Z(i,i) = sqrt(R::rchisq(df-i));
  }
  
  // Fill the lower matrix with random guesses
  for(unsigned int j = 0; j < m; j++){  
    for(unsigned int i = j+1; i < m; i++){    
      Z(i,j) = R::rnorm(0,1);
    }
  }
  
  // Lower triangle * chol decomp
  arma::mat C = arma::trimatl(Z).t() * arma::chol(S);
  
  // Return random wishart
  return C.t()*C;
}


//' @title Generate Random Inverse Wishart Distribution
//' @description Creates a random inverse wishart distribution when given degrees of freedom and a sigma matrix. 
//' @param df An \code{int} that represents the degrees of freedom.  (> 0)
//' @param S A \code{matrix} with dimensions m x m that provides Sigma, the covariance matrix. 
//' @return A \code{matrix} that is an inverse wishart distribution.
//' @seealso \code{\link{rwishart}}
//' @author James J Balamuta
//' @examples 
//' #Call with the following data:
//' riwishart(3, diag(2))
// [[Rcpp::export]]
arma::mat riwishart(unsigned int df, const arma::mat& S){
  return rwishart(df,S.i()).i();
}



// [[Rcpp::export]]

arma::mat GWishart_NOij_Gibbs_parallel( double bG, arma::mat DGr, arma::mat adjr, arma::mat Cr, double i, double j, double edgeij, int burnin, int nmc ) {
	
	//Function rWishart("rWishart");
	
	double p = DGr.n_rows;
	
	arma::mat DG (DGr.begin(),p,p);
//	NumericMatrix cppadj = clone(adjr);
	arma::mat adj (adjr.begin(),p,p);
//	NumericMatrix cppC = clone(Cr);
	arma::mat C (Cr.begin(),p,p);
	
	if( edgeij == 0 ){
		
		C(i,j) = 0;
		C(j,i) = 0;
		
		arma::mat DGjj;
		DGjj = DG(j,j);
		arma::mat l = arma::inv(DGjj);
		
		//arma::cube tempA = as<arma::cube>(rWishart(1,bG,wrap(l)));
		//arma::mat A = tempA.slice(0);
		arma::mat A = arma::wishrnd(l,bG);
		//arma::mat A = rwishart(bG,l);
		
		arma::mat C_12 = C.row(j);
		C_12.shed_col(j);
		
		arma::mat C_22 = C;
		C_22.shed_row(j);
		C_22.shed_col(j);
		arma::mat invC_22 = arma::inv(C_22);
		
		arma::mat c = C_12*invC_22*arma::trans(C_12);
		
		C(j,j) = A(0,0) + c(0,0);
		
		
	}else{
		
		arma::mat reorder = arma::linspace<arma::mat>(0, p-1, p);
		reorder.shed_row(i);
		reorder.shed_row(j-1);
		arma::mat ij = {i,j};
		reorder.insert_rows(reorder.size(),arma::trans(ij));
		
		arma::mat C_rows = C.row(reorder(0));
		for( int r = 1; r < p; r++ ){
			C_rows = arma::join_cols(C_rows,C.row(reorder(r)));
		}
		
		arma::mat C_reorder = C_rows.col(reorder(0));
		for( int c = 1; c < p; c++ ){
			C_reorder = arma::join_rows(C_reorder,C_rows.col(reorder(c)));
		}
		
		arma::mat R = arma::chol(C_reorder);
		double b = R(p-2,p-2);
		
		double m_post = -b*DG(i,j)/DG(j,j);
		double sig_post = 1/sqrt(DG(j,j));
		
		adj(i,j) = 1;
		adj(j,i) = 1;
		
		arma::vec rnorm = arma::randn<arma::vec>(1);
		R(p-2,p-1) = rnorm(0)*sig_post + m_post;
		
		double shape = bG/2;
		double scale = 2/DG(j,j);
		
		arma::vec rgamma = arma::randg<arma::vec>(1,arma::distr_param(shape,scale));
		R(p-1,p-1) = sqrt(rgamma(0));
		
		arma::mat C_updated = arma::trans(R.cols(p-2,p-1))*R.col(p-1);
		
		C(i,j) = C_updated(0);
		C(j,i) = C_updated(0);
		C(j,j) = C_updated(1);
		
	}
	
	arma::mat Sig = arma::inv(C);
	
	
	arma::uvec IsolatedNodeId_uvec  = find( arma::sum(adj,0) == 1 );
	arma::vec IsolatedNodeId  = arma::conv_to<arma::vec>::from(IsolatedNodeId_uvec);
	
	int INsize = IsolatedNodeId.n_elem;
	
	for( int iter = 0; iter < (burnin+nmc); iter++ ){
		
		if( INsize > 0 ){
			for( int i = 0; i < INsize; i++ ){
				
				double cliqueid = IsolatedNodeId(i);
				arma::mat DG_cliqueid;
				DG_cliqueid = DG(cliqueid,cliqueid);
				
				//arma::cube tempK_c = as<arma::cube>(rWishart(1,bG,wrap(arma::inv(DG_cliqueid))));
				//arma::mat K_c = tempK_c.slice(0);
				arma::mat K_c = arma::wishrnd(arma::inv(DG_cliqueid),bG);
				//arma::mat K_c = rwishart(bG,arma::inv(DG_cliqueid));
				C(cliqueid,cliqueid) = K_c(0,0);
				Sig(cliqueid,cliqueid) = 1/K_c(0,0);
				
			}
		}
		
		for( double i = 0; i < (p-1); i++ ){
			for( double j = i+1; j < p; j++ ){
				
				if( adj(i,j) == 1 ){
					
					double cliquesize = 2;
					
					arma::mat DG_row = arma::join_cols(DG.row(i),DG.row(j));
					arma::mat DG_cliqueid = arma::join_rows(DG_row.col(i),DG_row.col(j));
					arma::mat l = arma::inv(DG_cliqueid);
					l = (l+arma::trans(l))/2;
					arma::mat A = arma::wishrnd(l,bG+cliquesize-1);
					//arma::mat A = rwishart(bG+cliquesize-1,l);
					
					arma::mat C_12 = join_cols(C.row(i),C.row(j));
					C_12.shed_col(i);
					C_12.shed_col(j-1);
					
					arma::mat Sig_12 = join_cols(Sig.row(i),Sig.row(j));
					Sig_12.shed_col(i);
					Sig_12.shed_col(j-1);
					
					arma::mat Sig_22 = Sig;
					Sig_22.shed_row(i);
					Sig_22.shed_row(j-1);
					Sig_22.shed_col(i);
					Sig_22.shed_col(j-1);
					
					arma::mat Sig_row = join_cols(Sig.row(i),Sig.row(j));
					arma::mat Sig_11 = join_rows(Sig_row.col(i),Sig_row.col(j));
					arma::mat invSig_11 = arma::inv(Sig_11);
					invSig_11 = (invSig_11+arma::trans(invSig_11))/2;
					
					arma::mat invC_22 = Sig_22 - arma::trans(Sig_12)*invSig_11*Sig_12;
					
					arma::mat K_c = A + C_12*invC_22*arma::trans(C_12);
					K_c = (K_c + arma::trans(K_c))/2;
					
					arma::mat C_row = join_cols(C.row(i),C.row(j));
					arma::mat C_cliqueid = join_rows(C_row.col(i),C_row.col(j));
					
					arma::mat Delta = arma::inv(C_cliqueid-K_c);
					
					C(i,i) = K_c(0,0);
					C(i,j) = K_c(0,1);
					C(j,i) = K_c(1,0);
					C(j,j) = K_c(1,1);
					
					arma::mat aa = arma::inv(Delta-Sig_11);
					aa = ( aa + arma::trans(aa) )/2;
					
					arma::mat Sig_col = join_rows(Sig.col(i),Sig.col(j));
					arma::mat Sig_temp = Sig + Sig_col*aa*arma::trans(Sig_col);
					Sig = Sig_temp;
					
				}
				
			}
		}
		
	}
	
	
//	List Result;
//	Result["C"] = C;
//	Result["Sig"] = Sig;
	
	return(C);
	
}



// [[Rcpp::export]]

arma::mat GWishart_BIPS_pairwise_parallel( double bG, arma::mat DGr, arma::mat adjr, arma::mat Cr) {
	
	//Function rWishart("rWishart");
	
	double p = DGr.n_rows;
	
//	NumericMatrix cppD = clone(DGr);
	arma::mat DG (DGr.begin(),p,p);
	
//	NumericMatrix cppadj = clone(adjr);
	arma::mat adj (adjr.begin(),p,p);
	
//	NumericMatrix cppC = clone(Cr);
	arma::mat C (Cr.begin(),p,p);
	
	C = C%adj;
	arma::mat Sig = arma::inv(C);
	
	arma::uvec IsolatedNodeId_uvec  = find( arma::sum(adj,0) == 1 );
	arma::vec IsolatedNodeId  = arma::conv_to<arma::vec>::from(IsolatedNodeId_uvec);
	
	int INsize = IsolatedNodeId.n_elem;
		
	if( INsize > 0 ){
		for( int i = 0; i < INsize; i++ ){
			
			double cliqueid = IsolatedNodeId(i);
			
			arma::mat DG_cliqueid;
			DG_cliqueid = DG(cliqueid,cliqueid);
			arma::mat invDG_cliqueid = arma::inv(DG_cliqueid);
			
			//arma::cube tempK_c = as<arma::cube>(rWishart(1,bG,wrap(invDG_cliqueid)));
			//arma::mat K_c = tempK_c.slice(0);
			//arma::mat K_c = arma::wishrnd(invDG_cliqueid,bG);
			arma::mat K_c = rwishart(bG,invDG_cliqueid);
			
			C(cliqueid,cliqueid) = K_c(0,0);
			Sig(cliqueid,cliqueid) = 1/K_c(0,0);
			
		}
	}
	
	for( int i = 0; i < (p-1); i++ ){
		for( int j = (i+1); j < p; j++ ){
			
			if( adj(i,j) == 1 ){
				
				double cliquesize = 2;
				
				arma::mat DG_row = arma::join_cols(DG.row(i),DG.row(j));
				arma::mat DG_cliqueid = arma::join_rows(DG_row.col(i),DG_row.col(j));
				arma::mat invDG_cliqueid = arma::inv(DG_cliqueid);
				
		//		invDG_cliqueid = (invDG_cliqueid + arma::trans(invDG_cliqueid))/2;
				//arma::cube tempA = as<arma::cube>(rWishart(1,bG+cliquesize-1,wrap(invDG_cliqueid)));
				//arma::mat A = tempA.slice(0);
				//arma::mat A = arma::wishrnd(invDG_cliqueid,bG+cliquesize-1);
				arma::mat A = rwishart(bG+cliquesize-1,invDG_cliqueid);
				
				arma::mat C_12 = arma::join_cols(C.row(i),C.row(j));
				C_12.shed_col(i);
				C_12.shed_col(j-1);
				
				arma::mat Sig_12 = arma::join_cols(Sig.row(i),Sig.row(j));
				Sig_12.shed_col(i);
				Sig_12.shed_col(j-1);
				
				arma::mat Sig_22 = Sig;
				Sig_22.shed_row(i);
				Sig_22.shed_row(j-1);
				Sig_22.shed_col(i);
				Sig_22.shed_col(j-1);
				
				arma::mat Sig_row = arma::join_cols(Sig.row(i),Sig.row(j));
				arma::mat Sig_11 = join_rows(Sig_row.col(i),Sig_row.col(j));
				arma::mat invSig_11 = arma::inv(Sig_11);
				invSig_11 = (invSig_11 + arma::trans(invSig_11))/2;
				arma::mat invC_22 = Sig_22 - arma::trans(Sig_12)*invSig_11*Sig_12;
				
				arma::mat K_c = A + C_12*invC_22*arma::trans(C_12); 
				K_c = (K_c + arma::trans(K_c))/2;
				
				arma::mat C_row = arma::join_cols(C.row(i),C.row(j));
				arma::mat C_cliqueid = arma::join_rows(C_row.col(i),C_row.col(j));
				
				arma::mat Delta = arma::inv(C_cliqueid-K_c);
				
				C(i,i) = K_c(0,0);
				C(i,j) = K_c(0,1);
				C(j,i) = K_c(1,0);
				C(j,j) = K_c(1,1);
				
				arma::mat aa = arma::inv(Delta-Sig_11);
				aa = (aa + arma::trans(aa))/2;
				
				arma::mat Sig_col = arma::join_rows(Sig.col(i),Sig.col(j));
				
				Sig = Sig + Sig_col*aa*arma::trans(Sig_col);
				
			}
			
		}
	}
		
//	List Result;
//	Result["C"] = wrap(C);
//	Result["Sig"] = wrap(Sig);
	
	return(C);
	
}




// [[Rcpp::export]]

arma::cube GWishart_PAS_DMH_PG( double b_prior, arma::mat D_priorr, double n, arma::mat Sr, arma::mat Cr, arma::mat Betar, int burnin, int nmc ) {
	
	double p = D_priorr.n_rows;
	double nedge = p*(p-1)/2;
	
	arma::mat D_prior (D_priorr.begin(),p,p);
	arma::mat S (Sr.begin(),p,p);
	arma::mat Beta (Betar.begin(),p,p);
	
	double b_post = b_prior + n;
	arma::mat D_post = D_prior + S;
	
	arma::mat C (Cr.begin(),p,p);
	
//	arma::cube C_save(p,p,burnin+nmc+1);
//	arma::cube adj_save(p,p,burnin+nmc+1);
	
	arma::umat adj_umat = arma::abs(C)>(0.00005*arma::ones(p,p));
	arma::mat adj = arma::conv_to<arma::mat>::from(adj_umat);
	C = C%adj;
	
//	C_save.slice(0) = C;
//	adj_save.slice(0) = adj;
	
	arma::mat C_save (nedge+p,burnin+nmc+1);
	arma::mat adj_save (nedge+p,burnin+nmc+1);
	
	C_save.col(0) = upper2vec(C,1);
	adj_save.col(0) = upper2vec(adj,1);
	
	for( int iter = 0; iter < (burnin+nmc); iter++ ){
		
		for( double i = 0; i < (p-1); i++ ){
			for( double j = (i+1); j < p; j++ ){
				
				//Rcout << "i = " << i << std::endl;
				//Rcout << "j = " << i << std::endl;
				
				double w = log_H_parallel(b_prior,D_prior,n,S,(C),i,j)+ Beta(i,j);
				
				w = 1/(exp(w)+1);
				
				double current_ij = adj(i,j);
				
				arma::vec runif = arma::randu<arma::vec>(1);
				double propose_ij = runif(0)<w;
				
				if( propose_ij != current_ij ){
					
//					List temp = GWishart_NOij_Gibbs(b_prior,D_prior,wrap(adj),wrap(C),i,j,propose_ij,0,1);
//					arma::mat C_prop = as<arma::mat>(temp["C"]);
					arma::mat C_prop = GWishart_NOij_Gibbs_parallel(b_prior,D_prior,(adj),(C),i,j,propose_ij,0,1);
					
					if( det(C_prop) == 0.0 ){
						Rcout << "GWishart_NOij_Gibbs_parallel C_prop singular" << std::endl;
					}
					
					
					arma::vec eigval;
					arma::eig_sym(eigval,C_prop);
					if( eigval.min() <= 0.0 ){
						Rcout << "GWishart_NOij_Gibbs_parallel C_prop not PD" << std::endl;
					}
					
					double r0 = log_GWishart_NOij_pdf_parallel(b_prior,D_prior,(C_prop),i,j,current_ij);
					double r1 = log_GWishart_NOij_pdf_parallel(b_prior,D_prior,(C_prop),i,j,propose_ij);
					double r2 = r0 - r1;
					
					arma::vec runif1 = arma::randu<arma::vec>(1);
					arma::vec logrunif = log(runif1);
					
					if( logrunif(0) < r2 ){
						
						adj(i,j) = propose_ij;
						adj(j,i) = propose_ij;
						current_ij = propose_ij;
						
					}
					
				}
				
//				List temp1 = GWishart_NOij_Gibbs(b_post,D_post,wrap(adj),wrap(C),i,j,current_ij,0,0);
//				C = as<arma::mat>(temp1["C"]);
				arma::mat C_temp = GWishart_NOij_Gibbs_parallel(b_post,D_post,(adj),(C),i,j,current_ij,0,0);
				C = C_temp;
				
				if( det(C) == 0.0 ){
					Rcout << "GWishart_NOij_Gibbs_parallel C singular" << std::endl;
				}
				
				arma::vec eigval;
				arma::eig_sym(eigval,C);
				if( eigval.min() <= 0.0 ){
					Rcout << "GWishart_NOij_Gibbs_parallel C not PD" << std::endl;
				}
				
			}
		}
		
//		List temp2 = GWishart_BIPS_pairwise(b_post,wrap(D_post),wrap(adj),wrap(C));
//		C = as<arma::mat>(temp2["C"]);
		C = GWishart_BIPS_pairwise_parallel(b_post,(D_post),(adj),(C));
		
		if( det(C) == 0.0 ){
			Rcout << "GWishart_BIPS_pairwise_parallel C singular" << std::endl;
		}
		
		arma::vec eigval;
		arma::eig_sym(eigval,C);
		if( eigval.min() <= 0.0 ){
			Rcout << "GWishart_BIPS_pairwise_parallel C not PD" << std::endl;
		}
		
//		C_save.slice(iter) = C;
//		adj_save.slice(iter) = adj;
		
		C_save.col(iter+1) = upper2vec(C,1);
		adj_save.col(iter+1) = upper2vec(adj,1);
		
	}
	
//	List Result;
//	Result["C"] = wrap(C_save);
//	Result["adj"] = wrap(adj_save);
	
//	arma::cube Result = arma::join_slices(C_save,adj_save);
	
	arma::cube Result (nedge+p,burnin+nmc+1,2);
	Result.slice(0) = C_save;
	Result.slice(1) = adj_save;
	
	return(Result);
	
}




// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace RcppParallel;

struct GWishart_PAS_DMH_PG_worker : public Worker
{
	// source matrix
	double b_prior_iter;
	const RMatrix<double> D_prior_iter;
	double n_iter;
	
	int burninr;
	int nmcr;
	
	const RMatrix<double> data_p;
	const RMatrix<double> C_iterp;
	const RMatrix<double> Beta_iterp;
	
	// destination matrix
	RMatrix<double> C_save;
	RMatrix<double> adj_save;
	
	// initialize with source and destination
	GWishart_PAS_DMH_PG_worker( double b_prior_iter, const NumericMatrix& D_prior_iter, double n_iter, 
	const NumericMatrix& data_p, const NumericMatrix& C_iterp, const NumericMatrix& Beta_iterp, 
	int burninr, int nmcr, NumericMatrix C_save, NumericMatrix adj_save )
      : b_prior_iter(b_prior_iter), D_prior_iter(D_prior_iter), n_iter(n_iter), data_p(data_p), 
	  C_iterp(C_iterp), Beta_iterp(Beta_iterp), burninr(burninr), nmcr(nmcr), C_save(C_save), adj_save(adj_save) {}
	
	
	void operator()(std::size_t begin, std::size_t end) {
		
		//Function rWishart("rWishart");
		
		double b_prior = b_prior_iter;
		double n = n_iter;
		
		int burnin = burninr;
		int nmc = nmcr;
		
		double p = D_prior_iter.nrow();
		double numfiles = C_iterp.ncol();
		double nedge = p*(p-1)/2;
		
		RMatrix<double> cppD_prior = D_prior_iter;
		arma::mat D_prior (cppD_prior.begin(),p,p);
		
//		RMatrix<double> cppS_p = S_p;
//		arma::mat S_iter (cppS_p.begin(),nedge+p,numfiles,false);
		
		RMatrix<double> cppdata_p = data_p;
		arma::mat data (cppdata_p.begin(),n*numfiles,p);
		
		RMatrix<double> cppC_iter = C_iterp;
		arma::mat C_iter (cppC_iter.begin(),nedge+p,numfiles);
		
		RMatrix<double> cppBeta = Beta_iterp;
		arma::mat Beta_iter (cppBeta.begin(),nedge,numfiles);
		
		
		for( std::size_t l = begin; l < end; l++ ){
			
			arma::mat S (p,p,arma::fill::zeros);
//			arma::mat S_up (p,p,arma::fill::zeros);
//			arma::mat S_down (p,p,arma::fill::zeros);
			
//			S_up = vec2upper(S_iter.col(l),1);
//			S_down = arma::trans(S_up);
//			S_down.diag().zeros();
			
//			S = S_up + S_down;
			
			if( n == 1 ){
				S = (arma::trans(data.row(l))*data.row(l));
			}else{
				S = (arma::trans(data.rows(n*l,n*(l+1)-1))*data.rows(n*l,n*(l+1)-1));
			}
			
			arma::mat C (p,p,arma::fill::zeros);
			arma::mat C_up (p,p,arma::fill::zeros);
			arma::mat C_down (p,p,arma::fill::zeros);
			
			C_up = vec2upper(C_iter.col(l),1);
			C_down = arma::trans(vec2upper(C_iter.col(l),1));
			C_down.diag().zeros();
			
			C = C_up + C_down;
			
//			Rcout<< "file = " << l << std::endl;
//			Rcout<< "C = \n" << C << std::endl;
			
			arma::mat Beta (p,p,arma::fill::zeros);
			Beta = vec2upper(Beta_iter.col(l),0);
			
			arma::cube result = GWishart_PAS_DMH_PG(b_prior,(D_prior), n, (S), (C), (Beta), burnin, nmc );
			
//			arma::cube C_array = temp["C"];
			arma::mat C_array = result.slice(0);
			arma::mat C_result (nedge+p,1,arma::fill::zeros);
			if( burnin+nmc == 1 ){
				C_result = C_array.col(1);
			}else{
				C_result = mean(C_array.cols(burnin+1,burnin+nmc),1);
			}
			
			arma::mat adj_array = result.slice(1);
			arma::mat adj_result (nedge+p,1,arma::fill::zeros);
			if( burnin+nmc == 1 ){
				adj_result = adj_array.col(1);
			}else{
				arma::mat adj_temp = mean(adj_array.cols(burnin+1,burnin+nmc),1);
				arma::umat adj_umat = arma::abs(adj_temp)>(0.5*arma::ones(adj_temp.n_rows,adj_temp.n_cols));
				adj_result = arma::conv_to<arma::mat>::from(adj_umat);
			}
			
//			C_result = C_result%adj_result;
			
			for( int i = 0; i < (nedge+p); i++ ){
				double c = C_result(i);
				double a = adj_result(i);
				C_save(i,l) = c;
				adj_save(i,l) = a;
			}
			
		}
		
	}
	
};






// [[Rcpp::export]]

List GWishart_PAS_DMH_PG_parallel( double b_prior_iter, NumericMatrix D_prior_iter, 
	double n_iter, NumericMatrix data_p, NumericMatrix C_iterp, NumericMatrix Beta_iterp, 
	int burninr, int nmcr ){
	
	// allocate the output matrix
	NumericMatrix C_save(C_iterp.nrow(),C_iterp.ncol());	// C_iterp (nedge+p, numfiles)
	NumericMatrix adj_save(C_iterp.nrow(),C_iterp.ncol());
	
	// (pass input and output matrixes)
	GWishart_PAS_DMH_PG_worker GWishart_PAS_DMH_PG_w(b_prior_iter, D_prior_iter, n_iter,
	data_p, C_iterp, Beta_iterp, burninr, nmcr, C_save, adj_save);
	
	
	// call parallelFor to do the work
	parallelFor(0, C_iterp.ncol(), GWishart_PAS_DMH_PG_w);
	
	// return the output matrix
	List result;
	result["C"] = C_save;
	result["adj"] = adj_save;
	return result;
	
}






// [[Rcpp::export]]


double piecewise_coef( double x, double t, int n ){
	
	double a;
	
	if( x <= t ){
		a = arma::datum::pi*((n+1)/2)*(std::sqrt((2/(arma::datum::pi*x))*(2/(arma::datum::pi*x))*(2/(arma::datum::pi*x))))*exp(-(2*((n+1)/2)*((n+1)/2))/x);
	}else{
		a = arma::datum::pi*((n+1)/2)*exp(-((n+1)/2)*((n+1)/2)*((arma::datum::pi)*(arma::datum::pi))*x/2);
	}
	
	return a;
	
}



// [[Rcpp::export]]


double InverseGaussian_cdf( double x, double mu, double lambda ){
	
	double z1 = sqrt(lambda/x)*(x/mu-1);
	double z2 = -sqrt(lambda/x)*(x/mu+1);
	
	double cdf = arma::normcdf(z1)+exp(2*lambda/mu)*arma::normcdf(z2);
	
	return cdf;
	
}




// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace RcppParallel;

struct logit_worker : public Worker
{
	// source matrix
	const RMatrix<double> X;
	const RMatrix<double> beta_save_iter;
	
	// destination matrix
	RVector<double> w;
	
	// initialize with source and destination
	logit_worker(const NumericMatrix& X, const NumericMatrix& beta_save_iter, NumericVector w)
      : X(X), beta_save_iter(beta_save_iter), w(w) {}
	
	void operator()(std::size_t begin, std::size_t end) {
		
		double npar = X.ncol();
		double numfiles = X.nrow();
		RMatrix<double>  cppX = X;
		arma::mat armaX (cppX.begin(),numfiles,npar,false);
		RMatrix<double> cppbeta = beta_save_iter;
		arma::mat beta (cppbeta.begin(),npar,1,false);
		
		for( std::size_t i = begin; i < end; i++ ){
			
			arma::mat Xrow = armaX.row(i);
			
			arma::mat z_temp = arma::abs(Xrow*beta);
			double z = z_temp(0,0)/2;
			
			//// Sample J*(1,z) by accept-reject algorithm
			
			double K = (z*z)/2+(arma::datum::pi*arma::datum::pi)/8;
			
			double pz = 2*exp(-z)*InverseGaussian_cdf(0.64,1/z,1);
			double qz = (arma::datum::pi*exp(-K*0.64))/(2*K);
			
			arma::vec runif = arma::randu<arma::vec>(1);
			double u = runif(0);
			
			double cppdo = 1;
			while( cppdo == 1 ){
				
				double J;
				
				if( u > (pz/(pz+qz)) ){		// sample J' from Truncated Exponential
					
					arma::vec rgamma = arma::randg<arma::vec>(1,arma::distr_param(1,1));
					J = 0.64 + rgamma(0)/K;
					
				}else{		// sample J' from Truncated Inverse Gaussian
					
					double m = 1/z;
					
					if( m > 0.64 ){
						
						arma::vec rgamma = arma::randg<arma::vec>(1,arma::distr_param(1,1));
						arma::vec rgamma1 = arma::randg<arma::vec>(1,arma::distr_param(1,1));
						double E = rgamma(0);
						double E1 = rgamma1(0);
						
						while( (E*E) > (2*E1/0.64) ){
							rgamma = arma::randg<arma::vec>(1,arma::distr_param(1,1));
							rgamma1 = arma::randg<arma::vec>(1,arma::distr_param(1,1));
							E = rgamma(0);
							E1 = rgamma1(0);
						}
						
						J = 0.64/((1+0.64*E)*(1+0.64*E));
						double alpha = exp(-(z*z*J)/2);
						
						arma::vec runif1 = arma::randu<arma::vec>(1);
						double v = runif1(0);
						
						while( v > alpha ){
							rgamma = arma::randg<arma::vec>(1,arma::distr_param(1,1));
							rgamma1 = arma::randg<arma::vec>(1,arma::distr_param(1,1));
							E = rgamma(0);
							E1 = rgamma1(0);
							while( (E*E) > (2*E1/0.64) ){
								rgamma = arma::randg<arma::vec>(1,arma::distr_param(1,1));
								rgamma1 = arma::randg<arma::vec>(1,arma::distr_param(1,1));
								E = rgamma(0);
								E1 = rgamma1(0);
							}
							J = 0.64/((1+0.64*E)*(1+0.64*E));
							alpha = exp(-(z*z*J)/2);
							runif1 = arma::randu<arma::vec>(1);
							v = runif1(0);
						}
						
					}else{
						
						arma::vec runif2 = arma::randu<arma::vec>(1);
						double Y = runif2(0)*runif2(0);
						
						double J0 = m + 0.5*m*m*Y - 0.5*m*sqrt(4*m*Y + (m*Y)*(m*Y));
						
						arma::vec runif1 = arma::randu<arma::vec>(1);
						double v = runif1(0);
						
						if( v > m/(m+J0) ){
							J = m*m/J0;
						}else{
							J = J0;
						}
						
						while( J > 0.64 ){
							runif2 = arma::randu<arma::vec>(1);
							Y = runif2(0)*runif2(0);
							J0 = m + 0.5*m*m*Y - 0.5*m*sqrt(4*m*Y + (m*Y)*(m*Y));
							runif1 = arma::randu<arma::vec>(1);
							v = runif1(0);
							if( v > m/(m+J0) ){
								J = m*m/J0;
							}else{
								J = J0;
							}
						}
						
					}
					
				}
				
				int ind = 0;
				double S = piecewise_coef(J,0.64,ind);
				
				arma::vec runif1 = arma::randu<arma::vec>(1);
				double Y = S*runif1(0);
				
				while(1){
					ind = ind + 1;
					if( ind%2 == 1 ){
						S = S - piecewise_coef(J,0.64,ind);
						if( Y <= S ){
							w[i] = J/4;
							cppdo = 0;
							break;
						}
					}else{
						S = S + piecewise_coef(J,0.64,ind);
						if ( Y > S ){
							break;
						}
					}
					
				}
				
			}
		
		}
		
	}
	
	
};




// [[Rcpp::export]]

NumericVector logit_w_parallel(NumericMatrix X, NumericMatrix beta_save_iter){
	
	// allocate the output matrix
	NumericVector output(X.nrow());
	
	// (pass input and output matrixes)
	logit_worker logit_w(X, beta_save_iter, output);
	
	
	// call parallelFor to do the work
	parallelFor(0, X.nrow(), logit_w);
	
	// return the output matrix
	return output;
	
}



// [[Rcpp::export]]


arma::mat logit_parallel( NumericMatrix gijr, NumericMatrix Xr, int burnin, int nmc){
	
	double numfiles = Xr.nrow();
	double npar = Xr.ncol();
	NumericMatrix cppgij = clone(gijr);
	NumericMatrix cppX = clone(Xr);
	arma::mat gij (cppgij.begin(),numfiles,1,false);
	arma::mat X (cppX.begin(),numfiles,npar,false);
	
	//NumericMatrix cppmeanab = clone(meanabr);
	//arma::mat meanab (cppmeanab.begin(),npar,1,false);
	arma::mat meanab = arma::zeros(npar,1);
//	arma::mat meanab = arma::join_cols(arma::ones<arma::mat>(1,1)*log(0.15/0.85),arma::zeros(npar-1,1));
	
	arma::mat varab = arma::eye<arma::mat>(npar,npar)*100000000000000;
	arma::mat ivarab = arma::inv(varab);
	
	arma::mat beta_save = arma::zeros(npar,burnin+nmc+1);
	
	for( int iter = 0; iter < (burnin+nmc); iter++ ){
		
		
		//// Sample Polya-Gamma variable (-1z, with else, no const)
		
		arma::mat beta_col = beta_save.col(iter);
		
		NumericVector w_temp = logit_w_parallel(wrap(X),wrap(beta_col));
		arma::mat w (w_temp.begin(),numfiles,1,false);
		
		arma::mat W = arma::eye<arma::mat>(numfiles,numfiles);
		W.diag() = w;
		arma::mat kappa = gij-0.5;
		arma::mat varw = arma::inv(arma::trans(X)*W*X + ivarab);
		arma::mat meanw = varw*(arma::trans(X)*kappa + ivarab*meanab);
		
		arma::vec eigval;
		arma::mat L;
		arma::mat D = arma::zeros(npar,npar);
		arma::eig_sym( eigval, L, varw );
		D.diag() = eigval;
		
		beta_save.col(iter+1) = meanw + L*sqrt(D)*arma::mvnrnd(arma::zeros(npar,1),arma::eye<arma::mat>(npar,npar));
		
		
	}
	
	return(beta_save);
	
}



// [[Rcpp::export]]

double poisson_rng( double lambda ){
	
	// Poisson generator based upon the inversion by sequential search:
	
	double x = 0;
	double p = exp(-lambda);
	double s = p;
	
	arma::vec runif = arma::randu<arma::vec>(1);
	double u = runif(0);
	
	while( u > s ){
		x = x + 1;
		p = p*lambda/x;
		s = s + p;
	}
	
	return(x);
	
}



// [[Rcpp::export]]


double proposal_b( double b_from ){
	
	double q;
	
	if( b_from > 3.0 ){
		q = 1.0/2.0;
	}else if( b_from > 2.0 ){
		q = 1.0;
	}
	
	return q;
	
}


// [[Rcpp::export]]


double sample_b( double b_from){
	
	double b_prop = b_from;
	
	// degree of freedom > 2 in G-Wishart
	if( b_prop > 3.0 ){
		
		arma::vec runifvec = arma::randu<arma::vec>(1);
		double runif = runifvec(0);
		
		if( runif <= 1.0/2.0 ){
			--b_prop;
		}else{
			++b_prop;
		}
		
	}else if( b_prop > 2.0 ){
		
		++b_prop;
		
	}else{
		Rcout << "Warning: degree of freedom > 2 in G-Wishart" << std::endl;
	}
	
	return b_prop;
	
}



// [[Rcpp::export]]

double dNorm_rng( double mu, double sigma ){
	
	arma::vec v = arma::randn<arma::vec>(1);
//	double dx = mu + sigma*((floor(v(0))+0.5)/sqrt(1.0828));
	double dx = floor(mu+sigma*v(0));
	
	return(dx);
	
}




// [[Rcpp::export]]

double dNorm_pmf( double x, double mu, double sigma ){
	
//	double px = arma::normcdf(((x-mu)/sigma)*sqrt(1.0828)+0.5) - arma::normcdf(((x-mu)/sigma)*sqrt(1.0828)-0.5);
	double px = arma::normcdf((x+1-mu)/sigma) - arma::normcdf((x-mu)/sigma);
	
	return(px);
	
}



// [[Rcpp::export]]

arma::mat IPS( arma::mat L, arma::mat adj ){
	
	double p = L.n_rows;
	arma::mat K = arma::eye<arma::mat>(p,p);
	
	for( int i = 0; i < (p-1); i++ ){
		for( int j = (i+1); j < p; j++ ){
			
//			Rcout << "i =  " << i << std::endl;
//			Rcout << "j =  " << j << std::endl;
//			Rcout << "K before update\n" << K << std::endl;
			
			if( adj(i,j) == 1 ){
				arma::mat Ar = arma::join_cols(L.row(i),L.row(j));
				arma::mat A = arma::join_rows(Ar.col(i),Ar.col(j));
				
				arma::mat K11r = arma::join_cols(K.row(i),K.row(j));
				arma::mat K11 = arma::join_rows(K11r.col(i),K11r.col(j));
				
				arma::mat K12 = arma::join_cols(K.row(i),K.row(j));
				K12.shed_col(j);
				K12.shed_col(i);
				
				arma::mat K22 = K;
				K22.shed_row(j);
				K22.shed_row(i);
				K22.shed_col(j);
				K22.shed_col(i);
				
				arma::mat K21 = arma::trans(K12);
				
				arma::mat MK11 = arma::inv(A) + K12*arma::inv(K22)*K21;
				
				K(i,i) = MK11(0,0);
				K(i,j) = MK11(0,1);
				K(j,i) = MK11(1,0);
				K(j,j) = MK11(1,1);
				
			}
			
//			Rcout << "K after update\n" << K << std::endl;
			
		}
	}
	
	return(K);
	
}




// [[Rcpp::export]]


List BRUG( double b_priorr, NumericMatrix D_priorr, double n, NumericMatrix datar, NumericMatrix covariate, NumericMatrix Cr, 
NumericMatrix Betar, int burnin, int nmc, int init_burnin, int init_nmc, int inloop_burnin, int inloop_nmc, int PG_interval, int bD_interval ){
	
	double p = Cr.nrow();
	double nedge = p*(p-1)/2;
	double numfiles = datar.nrow();
	numfiles = numfiles/n;
	double npar = covariate.ncol();
	
	double b_prior = b_priorr;
	NumericMatrix cppD = clone(D_priorr);
	arma::mat D_prior (cppD.begin(),p,p,false);
	
//	D_prior = D_prior%arma::eye<arma::mat>(p,p);	// diagonalize D
	
	NumericMatrix cppC = clone(Cr);
	arma::mat C (cppC.begin(),p,p,false);
	NumericMatrix cppdata = clone(datar);
	arma::mat data (cppdata.begin(),numfiles*n,p,false);
	NumericMatrix cppX = clone(covariate);
	arma::mat X (cppX.begin(),numfiles,npar,false);
	NumericMatrix cppBeta = clone(Betar);
	arma::mat Beta (cppBeta.begin(),nedge,numfiles,false);
	
	arma::cube PG_beta (nedge,numfiles,burnin+nmc+1);
	arma::cube C_save (nedge+p,numfiles,burnin+nmc+1);	// inlcude_diag = 1
	arma::cube adj_save (nedge+p,numfiles,burnin+nmc+1);
	
	Rcout << "b = " << b_prior << std::endl;
	Rcout << "D = \n" << D_prior << std::endl;
	
	if( init_burnin+init_nmc > 0 ){
		
		List init = GWishart_PAS_DMH_PG_parallel(b_prior,wrap(D_prior), n, wrap(data), wrap(arma::repmat(upper2vec(C,1),1,numfiles)), wrap(Beta), init_burnin, init_nmc );
		C_save.slice(0) = as<arma::mat>(init["C"]);
		adj_save.slice(0) = as<arma::mat>(init["adj"]);
		
		
		arma::mat gij_init = adj_save.slice(0);
		for( double r = p; r > 0; r--){
			gij_init.shed_row(r*(r+1)/2-1);
		}
		
		arma::mat PG_beta_init (nedge,numfiles,arma::fill::zeros);
		arma::mat beta_init = arma::zeros(nedge,npar);
		
		for ( double e = 0; e < nedge; e++ ){
			
			arma::mat gijrow = gij_init.row(e);
			arma::mat beta_temp = logit_parallel(wrap(gijrow),wrap(X),500,1000);
			beta_init.row(e) = arma::trans(arma::mean(beta_temp,1));
			
			for( double l = 0; l < numfiles; l++ ){
				arma::mat pg = -beta_init.row(e)*arma::trans(X.row(l));
				PG_beta_init(e,l) = pg(0);
			}
			
		}
		
		PG_beta.slice(0) = PG_beta_init;
	}else{
		C_save.slice(0) = arma::repmat(upper2vec(C,1),1,numfiles);
		arma::umat adj_umat = (abs(C)>0);
		adj_save.slice(0) = arma::repmat(upper2vec(arma::conv_to<arma::mat>::from(adj_umat),1),1,numfiles);
		PG_beta.slice(0) = Beta;
	}
	
	arma::mat PG_beta_iter = PG_beta.slice(0);
	
	arma::vec b_prop (burnin+nmc+1);
	arma::cube D_prop (p,p,burnin+nmc+1);
	b_prop(0) = b_prior;
//	D_prop.slice(0) = D_prior;
	
	arma::mat sum_adj = arma::sum(adj_save.slice(0),1);
	arma::umat adj_common_umat = (sum_adj==(numfiles*arma::ones(sum_adj.n_rows,sum_adj.n_cols)));
	arma::mat adj_common = vec2upper(arma::conv_to<arma::mat>::from(adj_common_umat),1);
	adj_common = adj_common + arma::trans(adj_common) - arma::eye<arma::mat>(p,p);
	
	D_prop.slice(0) = arma::inv(IPS(D_prior,adj_common));
	
	double eta_D = p*500;	// df (>= p+1) parameter in Wishart proposal for D
	double eta = p+2;
	arma::mat V = arma::eye<arma::mat>(p,p)/eta;
	//double nu_D = p*480;	// df (> p+1) parameter in Inverse-Wishart proposal for D
	//double nu = p+2;
	//arma::mat Psi = arma::eye<arma::mat>(p,p);
	
	double lambda = p+4;
	
	//  sample b from posterior by MH
	
	double b0 = b_prop(0);
	arma::mat D0 = D_prop.slice(0);
	
	int accept_b = 0;
	int accept_D = 0;
	
	double b1 = sample_b(b0);
	
//	arma::cube C_k (p,p,numfiles);
//	arma::cube adj_k (p,p,numfiles);
	arma::mat C_iter = C_save.slice(0);
	arma::mat adj = adj_save.slice(0);
	
	
//	double logIG_b = 0;
	arma::mat sumC_j = arma::zeros(p,p);
	double logIG_b00 = 0;
	
	for( int k = 0; k < numfiles; k++ ){
		
//		arma::mat adj_temp = vec2upper(adj_save.slice(0).col(k),1);
//		adj_temp = adj_temp + arma::trans(adj_temp) - arma::eye<arma::mat>(p,p);
//		adj_k.slice(k) = adj_temp;
		
		//// log_hiwishart_InvA_const(makedecompgraph(adj),b,D) = -gwish_IG(adj,b,D) = -log(I_G(b,D))
		arma::mat C_up (p,p,arma::fill::zeros);
		arma::mat C_down (p,p,arma::fill::zeros);
		C_up = vec2upper(C_iter.col(k),1);
		C_down = arma::trans(C_up);
		C_down.diag().zeros();
		arma::mat C_temp = C_up + C_down;
//		C_k.slice(k) = C_temp;
		
//		logIG_b = logIG_b - log_hiwishart_InvA_const_parallel(adj_temp,b0,D0) + log_hiwishart_InvA_const_parallel(adj_temp,b1,D0) + (b1-b0)*log(arma::det(C_temp))/2;
//		logIG_b = logIG_b + log_GWishart_ud_const_mc_serial(b0,D0,adj_temp,100) - log_GWishart_ud_const_mc_serial(b1,D0,adj_temp,100) + (b1-b0)*log(arma::det(C_temp))/2;
		logIG_b00 = logIG_b00 + (b1-b0)*log(arma::det(C_temp))/2;
		sumC_j = sumC_j + C_temp;
		
	}
	
	arma::vec logIG_b0 = as<arma::vec>(log_GWishart_ud_const_mc_parallel(b0,wrap(D0),wrap(adj),100));
	arma::vec logIG_b1 = as<arma::vec>(log_GWishart_ud_const_mc_parallel(b1,wrap(D0),wrap(adj),100));
	double logIG_b = arma::accu(logIG_b0) - arma::accu(logIG_b1) + logIG_b00;
	
	double alpha_b = exp(logIG_b + lgamma(b0-3) - lgamma(b1-3) + (b1-b0)*log(lambda))*(proposal_b(b1)/proposal_b(b0));		// shifted pois prior, 2-way proposal
	
	arma::vec runif = arma::randu<arma::vec>(1);
	double vb = runif(0);
	
	if( vb <= alpha_b ){
		b_prop(0) = b1;
		b_prior = b1;
		++accept_b;
	}

	Rcout << "b = " << b_prior << std::endl;
	
	
	//  sample D from posterior by MH
	
	//arma::mat D1 = arma::wishrnd(D0/(eta_D),eta_D);
//	arma::mat D1 = rwishart(eta_D,D0/(eta_D));	// wishart proposal
	//arma::mat D1 = riwishart(nu_D,(nu_D-p-1)*D0);	// Inverse-wishart proposal
	
	arma::mat L = rwishart(eta_D,D0/(eta_D));	// wishart proposal
	arma::mat D1 = arma::inv(IPS(L,adj_common));
	
//	Rcout << "D1 = \n" << D1 << std::endl;
	
	
	arma::vec logIG_D0 = as<arma::vec>(log_GWishart_ud_const_mc_parallel(b_prior,wrap(D0),wrap(adj),100));
	arma::vec logIG_D1 = as<arma::vec>(log_GWishart_ud_const_mc_parallel(b_prior,wrap(D1),wrap(adj),100));
	double logIG_D = arma::accu(logIG_D0) - arma::accu(logIG_D1);
	
//	double constD = sqrt(std::pow(arma::det(D0)/arma::det(D1),2*eta_D-eta));	// wishart prior, wishart proposal
	double logconstD = (2*eta_D-eta)*(log(arma::det(D0))-log(arma::det(D0)))/2;	// wishart prior, wishart proposal
	//double constD = sqrt(std::pow(arma::det(D1)/arma::det(D0),2*nu_D-nu));	// iwishart prior, iwishart proposal
	//double constD = sqrt(std::pow(arma::det(D0)/arma::det(D1),2*eta_D+nu));	// iwishart prior, wishart proposal
	
	double traceD = arma::trace( (sumC_j + arma::inv(V))*(D1-D0) + arma::inv(D1/eta_D)*D0 - arma::inv(D0/eta_D)*D1 )/2;
	
	double alpha_D = exp(logIG_D + logconstD - traceD);	// wishart prior, wishart proposal
	
	arma::vec runif1 = arma::randu<arma::vec>(1);
	double vD = runif1(0);
	
	if( vD <= alpha_D ){
//		D_prop.slice(iter+1) = IPS(D1,adj.slice(0));
//		D_prop.slice(iter+1) = arma::inv(IPS(D1,adj.slice(0)));
		D_prop.slice(0) = D1;
		D_prior = D1;
		++accept_D;
	}else{
		D_prop.slice(0) = D0;
	}
	
	Rcout << "D = \n" << D_prior << std::endl;
	
//	Rcout << "The 1st graph is \n" << C_k.slice(0) << std::endl;
	
	
	
	for( int iter = 0; iter < (burnin+nmc); iter++ ){
		
		List result = GWishart_PAS_DMH_PG_parallel(b_prior,wrap(D_prior), n, wrap(data), wrap(C_save.slice(iter)), wrap(PG_beta_iter), inloop_burnin, inloop_nmc );
		
		C_iter = as<arma::mat>(result["C"]);
		C_save.slice(iter+1) = C_iter;
		adj = as<arma::mat>(result["adj"]);
		adj_save.slice(iter+1) = adj;
		
		
		
	//	if( iter - floor(iter/150.0)*150.0 == 0 ){
			
	//		Rcout << "*************** Iter = " << iter << std::endl;
	//		Rcout << "b = " << b_prior << std::endl;
	//		Rcout << "D = \n" << D_prior << std::endl;
			
			
	//		b0 = b_prop(iter);
			
	//		b1 = sample_b(b0);
			
	//		logIG_b = 0;
			
	//		for( int k = 0; k < numfiles; k++ ){
				//// log_hiwishart_InvA_const(makedecompgraph(adj),b,D) = -gwish_IG(adj,b,D) = -log(I_G(b,D))
			//	logIG_b = logIG_b - log_hiwishart_InvA_const_parallel(adj_k.slice(k),b0,D_prior) + log_hiwishart_InvA_const_parallel(adj_k.slice(k),b1,D_prior) + (b1-b0)*log(arma::det(C_k.slice(k)))/2;
	//			logIG_b = logIG_b + log_GWishart_ud_const_mc_serial(b0,D_prior,adj_k.slice(k),100) - log_GWishart_ud_const_mc_serial(b1,D_prior,adj_k.slice(k),100) + (b1-b0)*log(arma::det(C_k.slice(k)))/2;
	//		}
			
	//		alpha_b = exp(logIG_b + lgamma(b0-3) - lgamma(b1-3) + (b1-b0)*log(lambda))*(proposal_b(b1)/proposal_b(b0));		// shifted pois prior, 2-way proposal
			
	//		runif = arma::randu<arma::vec>(1);
	//		vb = runif(0);
			
	//		if( vb <= alpha_b ){
	//			b_prop(iter+1) = b1;
	//			b_prior = b1;
	//			++accept_b;
	//		}else{
	//			b_prop(iter+1) = b0;
	//			b_prior = b0;
	//		}
			
	//	}
		
		
		if( iter - floor(iter/bD_interval)*bD_interval == 0 ){
			
			
			sum_adj = arma::sum(adj,1);
			adj_common_umat = (sum_adj==(numfiles*arma::ones(sum_adj.n_rows,sum_adj.n_cols)));
			adj_common = vec2upper(arma::conv_to<arma::mat>::from(adj_common_umat),1);
			adj_common = adj_common + arma::trans(adj_common) - arma::eye<arma::mat>(p,p);
			
			//  sample b from posterior by MH
			
			b0 = b_prop(iter);
			D0 = arma::inv(IPS(D_prop.slice(iter),adj_common));
			
			b1 = sample_b(b0);
			
			sumC_j = arma::zeros(p,p);
			logIG_b00 = 0;
			
			for( int k = 0; k < numfiles; k++ ){
				
		//		arma::mat adj_temp = vec2upper(adj_save.slice(0).col(k),1);
		//		adj_temp = adj_temp + arma::trans(adj_temp) - arma::eye<arma::mat>(p,p);
		//		adj_k.slice(k) = adj_temp;
				
				//// log_hiwishart_InvA_const(makedecompgraph(adj),b,D) = -gwish_IG(adj,b,D) = -log(I_G(b,D))
				arma::mat C_up (p,p,arma::fill::zeros);
				arma::mat C_down (p,p,arma::fill::zeros);
				C_up = vec2upper(C_iter.col(k),1);
				C_down = arma::trans(C_up);
				C_down.diag().zeros();
				arma::mat C_temp = C_up + C_down;
		//		C_k.slice(k) = C_temp;
				
		//		logIG_b = logIG_b - log_hiwishart_InvA_const_parallel(adj_temp,b0,D0) + log_hiwishart_InvA_const_parallel(adj_temp,b1,D0) + (b1-b0)*log(arma::det(C_temp))/2;
		//		logIG_b = logIG_b + log_GWishart_ud_const_mc_serial(b0,D0,adj_temp,100) - log_GWishart_ud_const_mc_serial(b1,D0,adj_temp,100) + (b1-b0)*log(arma::det(C_temp))/2;
				logIG_b00 = logIG_b00 + (b1-b0)*log(arma::det(C_temp))/2;
				sumC_j = sumC_j + C_temp;
				
			}
			
			logIG_b0 = as<arma::vec>(log_GWishart_ud_const_mc_parallel(b0,wrap(D0),wrap(adj),100));
			logIG_b1 = as<arma::vec>(log_GWishart_ud_const_mc_parallel(b1,wrap(D0),wrap(adj),100));
			logIG_b = arma::accu(logIG_b0) - arma::accu(logIG_b1) + logIG_b00;
		
			alpha_b = exp(logIG_b + lgamma(b0-3) - lgamma(b1-3) + (b1-b0)*log(lambda))*(proposal_b(b1)/proposal_b(b0));		// shifted pois prior, 2-way proposal
			
			runif = arma::randu<arma::vec>(1);
			vb = runif(0);
			
			if( vb <= alpha_b ){
		//		b_prop(iter+1) = b1;
				b_prior = b1;
				++accept_b;
			}else{
		//		b_prop(iter+1) = b0;
				b_prior = b0;
			}
			
			
			//  sample D from posterior by MH
			
			//arma::mat D1 = arma::wishrnd(D0/(eta_D),eta_D);
		//	arma::mat D1 = rwishart(eta_D,D0/(eta_D));	// wishart proposal
			//arma::mat D1 = riwishart(nu_D,(nu_D-p-1)*D0);	// Inverse-wishart proposal
			
			L = rwishart(eta_D,D0/(eta_D));	// wishart proposal
			D1 = arma::inv(IPS(L,adj_common));
			
		//	Rcout << "D1 = \n" << D1 << std::endl;
			
			logIG_D0 = as<arma::vec>(log_GWishart_ud_const_mc_parallel(b_prior,wrap(D0),wrap(adj),100));
			logIG_D1 = as<arma::vec>(log_GWishart_ud_const_mc_parallel(b_prior,wrap(D1),wrap(adj),100));
			logIG_D = arma::accu(logIG_D0) - arma::accu(logIG_D1);
			
		//	constD = sqrt(std::pow(arma::det(D0)/arma::det(D1),2*eta_D-eta));	// wishart prior, wishart proposal
			logconstD = (2*eta_D-eta)*(log(arma::det(D0))-log(arma::det(D0)))/2;	// wishart prior, wishart proposal
			//constD = sqrt(std::pow(arma::det(D1)/arma::det(D0),2*nu_D-nu));	// iwishart prior, iwishart proposal
			//constD = sqrt(std::pow(arma::det(D0)/arma::det(D1),2*eta_D+nu));	// iwishart prior, wishart proposal
			
			traceD = arma::trace( (sumC_j + arma::inv(V))*(D1-D0) + arma::inv(D1/eta_D)*D0 - arma::inv(D0/eta_D)*D1 )/2;
			
			alpha_D = exp(logIG_D + logconstD - traceD);	// wishart prior, wishart proposal
			
			runif1 = arma::randu<arma::vec>(1);
			vD = runif1(0);
			
			if( vD <= alpha_D ){
		//		D_prop.slice(iter+1) = IPS(D1,adj.slice(0));
		//		D_prop.slice(iter+1) = arma::inv(IPS(D1,adj.slice(0)));
		//		D_prop.slice(iter+1) = D1;
				D_prior = D1;
				++accept_D;
			}else{
		//		D_prop.slice(iter+1) = D0;
				D_prior = D0;
			}
			
			
		//	Rcout << "*************** Iter = " << iter << std::endl;
		//	Rcout << "b = " << b_prior << std::endl;
		//	Rcout << "D = \n" << D_prior << std::endl;
		//	Rcout << "alpha_b = " << alpha_b << std::endl;
		//	Rcout << "alpha_D = " << alpha_D << std::endl;
			
		}
		
		if( iter - floor(iter/100)*100 == 0 ){
			Rcout << "*************** Iter = " << iter << std::endl;
			Rcout << "b = " << b_prior << std::endl;
			Rcout << "D = \n" << D_prior << std::endl;
		}
		
		if( iter - floor(iter/PG_interval)*PG_interval == 0 ){
			
			arma::mat gij;
			arma::mat beta_est = arma::zeros(nedge,npar);
			
			if( (iter/PG_interval) > 0 ){
				
		//		arma::mat adj_temp0 = mean(adj_save.slices(PG_interval*(iter/PG_interval-1)+2,iter+1),2);
		//		arma::umat adj_temp = adj_temp0>(0.5*arma::ones(adj_temp0.n_rows,adj_temp0.n_cols));
				
		//		gij = arma::conv_to<arma::mat>::from(adj_temp);
				gij = adj_save.slice(iter+1);
				
				for( double r = p; r > 0; r--){
					gij.shed_row(r*(r+1)/2-1);
				}
				
				for ( double e = 0; e < nedge; e++ ){
					
					arma::mat gijrow = gij.row(e);
					arma::mat beta_temp = logit_parallel(wrap(gijrow),wrap(X),500,1000);
					beta_est.row(e) = arma::trans(arma::mean(beta_temp,1));
					
					for( double l = 0; l < numfiles; l++ ){
						arma::mat pg = -beta_est.row(e)*arma::trans(X.row(l));
						PG_beta_iter(e,l) = pg(0,0);
					}
					
				}
				
			}
			
		}
		
		b_prop(iter+1) = b_prior;
		D_prop.slice(iter+1) = D_prior;		
		PG_beta.slice(iter+1) = PG_beta_iter;
		
	}
	
	List Result;
	Result["C"] = wrap(C_save);
	Result["adj"] = wrap(adj_save);
	Result["b"] = wrap(b_prop);
	Result["D"] = wrap(D_prop);
	Result["accept_b"] = accept_b;
	Result["accept_D"] = accept_D;
	
	
	return (Result);
	
}





















































