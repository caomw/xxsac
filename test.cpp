#include <theia/theia.h>

#include <Eigen/Dense>

#include <vector>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>

using namespace theia;
using namespace std;

// PointPair is the data structure used in estimation
struct PointPair{
	Eigen::Vector2d p1, p2;
	double scale1, scale2;
	double orientation1, orientation2;
	double score;
};

struct FundModel{
	Eigen::Matrix3d F;
};

class FundEstimator: public Estimator<PointPair, FundModel> {
	public:
	// number of point pairs needed to estimate the fundamental matrix
	double SampleSize() const { return 8; }

	// estimate the fundamental matrix from eight point-pairs
	bool EstimateModel(const vector<PointPair> &data, vector<FundModel> *models) const{
		FundModel model;
		vector< Eigen::Vector2d > kp1;
		vector< Eigen::Vector2d > kp2;
		for (int i = 0; i < 8; i++){
			kp1.push_back( data[i].p1 );
			kp2.push_back( data[i].p2 );
		}
		bool estimation_success = NormalizedEightPoint( kp1, kp2, &(model.F) );
		
		models->push_back( model );
		return estimation_success;
	}

	double Error(const PointPair &data, const FundModel &model) const {
		Eigen::Vector3d x1, x2;
		x1[0] = data.p1[0];
		x1[1] = data.p1[1];
		x2[0] = data.p2[0];
		x2[1] = data.p2[1];
		x1[2] = x2[2] = 1;

		double x2tFx1 = x2.transpose()*model.F*x1;
		Eigen::Vector3d Fx1 = model.F*x1;
		Eigen::Vector3d Ftx2 = model.F.transpose()*x2;
		double distance = x2tFx1*x2tFx1 / (Fx1[0]*Fx1[0] + Fx1[1]*Fx1[1] + Ftx2[0]*Ftx2[0] + Ftx2[1]*Ftx2[1]);
		return distance;
	}
};

bool get_putative_match_data( string file_name, vector< PointPair > &match ){
	ifstream kp_file( file_name );
	double x, y, scale, orientation;
	int num;
	vector< Eigen::Vector2d > kp1;
	vector< Eigen::Vector2d > kp2;
	if (kp_file.is_open()){
		kp_file >> num;
		if ( !match.empty() ){
			match.clear();
		}
		for (int i = 0; i < num; i++){
			kp_file >> x >> y >> scale >> orientation;
			kp1.push_back( Eigen::Vector2d(x, y) );
			match.push_back( PointPair() );
			match[i].p1 << x, y;
			match[i].scale1 = scale;
			match[i].orientation1 = orientation;
		}
		for (int i = 0; i < num; i++){
			kp_file >> x >> y >> scale >> orientation;
			kp2.push_back( Eigen::Vector2d(x, y) );
			match[i].p2 << x, y;
			match[i].scale2 = scale;
			match[i].orientation2 = orientation;
		}
		kp_file.close();
		return true;
	}else{
		return false;
	}
}

bool output_result(string filename, RansacSummary &summary, FundModel &model){
	ofstream outfile(filename);
        if (outfile.is_open()){
                outfile << model.F << endl;
		// iteration num
		outfile << summary.num_iterations << endl;
                outfile << summary.inliers.size() << endl;
                for (int i = 0; i < summary.inliers.size(); i++){
              	  outfile << summary.inliers.at(i) << " ";
                }
                outfile.close();
		return true;
        }else{
                cout << "error opening file" << endl;
		return false;
        }
}

bool get_F_matrix(string filename, Eigen::Matrix3d &F){
	ifstream infile( filename );
	if (infile.is_open()){
		for (int i = 0; i < 3; i++){
			for (int j = 0; j < 3; j++){
				infile >> F(i, j);
			}
		}
		infile.close();
		return true;
	}else{
		cout << "failed to open file " << filename << endl;
		return false;
	}
	cout << "test Fundamental matrix: " << endl;
	cout << F << endl;
}
void check_inliers(double Threshold, Eigen::Matrix3d &F, vector< PointPair > &match, RansacSummary &summary, vector< int > &trueInliers, int &trueInlierNum, double &inlierRatio){
	if ( !trueInliers.empty() ){
		trueInliers.clear();
	}
	trueInlierNum = 0;
	inlierRatio = 0;
	for (int i = 0; i < summary.inliers.size(); i++){
		//Eigen::Vector3d x1, x2;
		//x1[0] = match[ summary.inliers[i] ].p1[0];
		//x1[1] = match[ summary.inliers[i] ].p1[1];
		//x2[0] = match[ summary.inliers[i] ].p2[0];
		//x2[1] = match[ summary.inliers[i] ].p2[1];
		//x1[2] = x2[2] = 1;

		//double x2tFx1 = x2.transpose()*F*x1;
		//Eigen::Vector3d Fx1 = F*x1;
		//Eigen::Vector3d Ftx2 = F.transpose()*x2;
		//double distance = x2tFx1*x2tFx1 / (Fx1[0]*Fx1[0] + Fx1[1]*Fx1[1] + Ftx2[0]*Ftx2[0] + Ftx2[1]*Ftx2[1]);
		//cout << "test distance: " << distance << endl;
		//if (distance < Threshold){
		//	trueInliers.push_back(1);
		//	++trueInlierNum;
		//}else{
		//	trueInliers.push_back(0);
		//}
		Eigen::Vector2d x1( match[ summary.inliers[i] ].p1 );
		Eigen::Vector2d x2( match[ summary.inliers[i] ].p2 );
		Eigen::Vector3d x1h( x1[0], x1[1], 1);
		x1h = F*x1h;
		Eigen::Vector2d x1_new( x1h[0] / x1h[2], x1h[1] / x1h[2] );
		Eigen::Vector2d r(x2 - x1_new);
		if (sqrt(r[0]*r[0]+r[1]*r[1]) < Threshold){
			trueInliers.push_back(1);
			++trueInlierNum;
		}else{
			trueInliers.push_back(0);
		}
	}
	inlierRatio = trueInlierNum / (double)summary.inliers.size();
}
int main(){
	// define Fundamental matrix and Keypoints
	Eigen::Matrix3d F;
	vector< PointPair > match;


	double thresh = 2;
	double inlierRatio;
	int iterationNum;
	int trueInlierNum;
	int trialNum = 20;

	double totTrueInlierNum;
	double totInlierRatio;
	double totIterationNum;
	double totInlierNum;
	vector< int > trueInliers;


	// set the fundmatrix estimator
        FundEstimator fund_estimator;
	// basic ransac parameters
        RansacParameters params;
        params.error_thresh = 1;
	params.max_iterations = 10000;

	// read data from txt file
	string setname, filename;

	cout << "input dataset name" << endl;
	cin >> setname;
	cout << "=====" << endl;

	cout << "input F filename" << endl;
	cin >> filename;
	get_F_matrix("/media/sf_workspace/"+setname+"/"+filename, F);

	cout << "input summary filename" << endl;
	cin >> filename;
	ofstream summary_file( filename );
	if ( !summary_file.is_open() ){
		cout << "failed to open summary file" << endl;
		return -1;
	}
	


	//////////
	//PRUNE+RANSAC
	/////////
	cout << "Ransac putative PRUNED match file:" << endl;
	cin >> filename;
	get_putative_match_data("/media/sf_workspace/"+setname+"/"+filename, match);

	cout << "Ransac NON-pruned match file:" << endl;
	string filename2;
	cin >> filename2;
	vector< PointPair > nonPruneMatch;
	get_putative_match_data("/media/sf_workspace/"+setname+"/"+filename2, nonPruneMatch);



				int putativeInlier = 0;
				for (int j = 0; j < nonPruneMatch.size(); j++){
					Eigen::Vector2d x1( nonPruneMatch[j].p1 );
					Eigen::Vector2d x2( nonPruneMatch[j].p2 );
					Eigen::Vector3d x1h( x1[0], x1[1], 1);
					x1h = F*x1h;
					Eigen::Vector2d x1_new( x1h[0] / x1h[2], x1h[1] / x1h[2] );
					Eigen::Vector2d r(x2 - x1_new);
					if (sqrt(r[0]*r[0]+r[1]*r[1]) < thresh){
						++putativeInlier;
					}
				}
				cout << "Putative ture inlier ratio:" << putativeInlier / (double) nonPruneMatch.size() << endl;

	vector< int > ourInlier;
		

        cout << "Prune+Ransac running" << endl;
        Ransac<FundEstimator> prune_estimator(params, fund_estimator);

	FundModel prune_model;
        RansacSummary prune_summary;


	totTrueInlierNum = 0;
	totInlierRatio   = 0;
	totIterationNum  = 0;
	totInlierNum = 0;
	for (int i = 0; i < trialNum; i++){
		cout << i << endl;
		prune_estimator.Initialize();
		cout << "before" << endl;
		prune_estimator.Estimate(match, &prune_model, &prune_summary);
		cout << "end" << endl;
		//check_inliers(thresh, F, match, prune_summary, trueInliers, trueInlierNum, inlierRatio);
		//totTrueInlierNum += trueInlierNum;
		//totInlierRatio   += inlierRatio;

		int inlier_num = 0;
		int trueInlierNum = 0;
		if (!ourInlier.empty()){
			ourInlier.clear();
		}
		for (int j = 0; j < nonPruneMatch.size(); j++){
			if ( fund_estimator.Error( nonPruneMatch[j], prune_model ) < 1 ){
				++inlier_num;
				ourInlier.push_back(j);
				// check if this is a true inlier
				Eigen::Vector2d x1( nonPruneMatch[j].p1 );
				Eigen::Vector2d x2( nonPruneMatch[j].p2 );
				Eigen::Vector3d x1h( x1[0], x1[1], 1);
				x1h = F*x1h;
				Eigen::Vector2d x1_new( x1h[0] / x1h[2], x1h[1] / x1h[2] );
				Eigen::Vector2d r(x2 - x1_new);
				if (sqrt(r[0]*r[0]+r[1]*r[1]) < thresh){
					++trueInlierNum;
				}
			}
		}
		totTrueInlierNum += trueInlierNum;
		totInlierNum += inlier_num;
		totInlierRatio +=trueInlierNum/(double)inlier_num;
		totIterationNum  += prune_summary.num_iterations;
	}
	totIterationNum /= trialNum;
	totInlierRatio  /= trialNum;
	totTrueInlierNum/= trialNum;
	totInlierNum /= trialNum;
        cout << "iteration num:" << totIterationNum << endl;
	cout << "ture inlier num:" << totTrueInlierNum << endl;
	cout << "inlier num:" << totInlierNum << endl;
	cout << "inlier ratio:" << totInlierRatio << endl;
	ofstream outfile("/media/sf_workspace/"+setname+"/ransac_"+filename);
        if (outfile.is_open()){
		// iteration num
                for (int i = 0; i < ourInlier.size(); i++){
              	  outfile << ourInlier[i] << " ";
                }
                outfile.close();
        }else{
                cout << "error opening file" << endl;
        }
	cout << "=====" << endl;

	summary_file << "Prune Ransac" << endl;
        summary_file << "iteration num:" << totIterationNum << endl;
	summary_file << "ture inlier num:" << totTrueInlierNum << endl;
	summary_file << "inlier ratio:" << totInlierRatio << endl;
	summary_file << "=====" << endl << endl;
	cout << "ok" << endl;


	//////////
	//RANSAC
	/////////
	cout << "Ransac putative match file:" << endl;
	cin >> filename;
	get_putative_match_data("/media/sf_workspace/"+setname+"/"+filename, match);
		

        cout << "Ransac running" << endl;
        Ransac<FundEstimator> ransac_estimator(params, fund_estimator);

	FundModel ransac_model;
        RansacSummary ransac_summary;

	totTrueInlierNum = 0;
	totInlierRatio   = 0;
	totIterationNum  = 0;
	totInlierNum = 0;
	for (int i = 0; i < trialNum; i++){
		ransac_estimator.Initialize();
		ransac_estimator.Estimate(match, &ransac_model, &ransac_summary);
		//cout << ransac_summary.num_iterations << endl;
	check_inliers(thresh, F, match, ransac_summary, trueInliers, trueInlierNum, inlierRatio);
	totTrueInlierNum += trueInlierNum;
	totInlierRatio   += inlierRatio;
	totIterationNum  += ransac_summary.num_iterations;
	totInlierNum     += ransac_summary.inliers.size();
}
totIterationNum /= trialNum;
totInlierRatio  /= trialNum;
totTrueInlierNum/= trialNum;
totInlierNum    /= trialNum;
cout << "iteration num:" << totIterationNum << endl;
cout << "ture inlier num:" << totTrueInlierNum << endl;
cout << "inlier num:" << totInlierNum << endl;
cout << "inlier ratio:" << totInlierRatio << endl;
output_result("/media/sf_workspace/"+setname+"/ransac_"+filename, ransac_summary, ransac_model);
cout << "=====" << endl;

summary_file << " Ransac" << endl;
summary_file << "iteration num:" << totIterationNum << endl;
summary_file << "ture inlier num:" << totTrueInlierNum << endl;
summary_file << "inlier ratio:" << totInlierRatio << endl;
summary_file << "=====" << endl << endl;



/////////
//PROSAC
////////
cout << "Prosac putative match file:" << endl;
cin >> filename;
get_putative_match_data("/media/sf_workspace/"+setname+"/"+filename, match);

cout << "Prosac running" << endl;
Prosac<FundEstimator> prosac_estimator(params, fund_estimator);
prosac_estimator.Initialize();

RansacSummary prosac_summary;
FundModel prosac_model;

totTrueInlierNum = 0;
totInlierRatio   = 0;
totIterationNum  = 0;
totInlierNum = 0;
for (int i = 0; i < trialNum; i++){
	prosac_estimator.Initialize();
	prosac_estimator.Estimate(match, &prosac_model, &prosac_summary);
	//cout << prosac_summary.num_iterations << endl;
	check_inliers(thresh, F, match, prosac_summary, trueInliers, trueInlierNum, inlierRatio);
		totTrueInlierNum += trueInlierNum;
		totInlierRatio   += inlierRatio;
		totIterationNum  += prosac_summary.num_iterations;
		totInlierNum     += prosac_summary.inliers.size();
	}
	totIterationNum /= trialNum;
	totInlierRatio  /= trialNum;
	totTrueInlierNum/= trialNum;
	totInlierNum    /= trialNum;
        cout << "iteration num:" << totIterationNum << endl;
	cout << "ture inlier num:" << totTrueInlierNum << endl;
	cout << "inlier num:" << totInlierNum << endl;
	cout << "inlier ratio:" << totInlierRatio << endl;
	output_result("/media/sf_workspace/"+setname+"/prosac_"+filename, prosac_summary, prosac_model);
	cout << "=====" << endl;
    
	summary_file << " Prosac" << endl;
        summary_file << "iteration num:" << totIterationNum << endl;
	summary_file << "ture inlier num:" << totTrueInlierNum << endl;
	summary_file << "inlier ratio:" << totInlierRatio << endl;
	summary_file << "=====" << endl << endl;
        ////////
	//MLESAC
	////////
        cout << "Mlesac putative match file:" << endl;
   	cin >> filename; 
	get_putative_match_data("/media/sf_workspace/"+setname+"/"+filename, match);

	
	cout << "Mlesac running" << endl;
        Mlesac<FundEstimator> mlesac_estimator(params, fund_estimator);

        FundModel mlesac_model;
        RansacSummary mlesac_summary;

	totTrueInlierNum = 0;
	totInlierRatio   = 0;
	totIterationNum  = 0;
	totInlierNum     = 0;
	for (int i = 0; i < trialNum; i++){
		mlesac_estimator.Initialize();
		mlesac_estimator.Estimate(match, &mlesac_model, &mlesac_summary);
		//cout << mlesac_summary.num_iterations << endl;
		check_inliers(thresh, F, match, mlesac_summary, trueInliers, trueInlierNum, inlierRatio);
		totTrueInlierNum += trueInlierNum;
		totInlierRatio   += inlierRatio;
		totIterationNum  += mlesac_summary.num_iterations;
		totInlierNum     += mlesac_summary.inliers.size();
	}
	totIterationNum /= trialNum;
	totInlierRatio  /= trialNum;
	totTrueInlierNum/= trialNum;
	totInlierNum    /= trialNum;
        cout << "iteration num:" << totIterationNum << endl;
	cout << "ture inlier num:" << totTrueInlierNum << endl;
	cout << "inlier num:" << totInlierNum << endl;
	cout << "inlier ratio:" << totInlierRatio << endl;
	output_result("/media/sf_workspace/"+setname+"/mlesac_"+filename, mlesac_summary, mlesac_model);
	cout << "=====" << endl;

	summary_file << " Mlesac" << endl;
        summary_file << "iteration num:" << totIterationNum << endl;
	summary_file << "ture inlier num:" << totTrueInlierNum << endl;
	summary_file << "inlier ratio:" << totInlierRatio << endl;
	summary_file << "=====" << endl << endl;

	summary_file.close();

    

}

