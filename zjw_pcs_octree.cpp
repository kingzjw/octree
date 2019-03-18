#pragma once
#include "zjw_pcs_octree.h"

PcsOctree::PcsOctree(int _m, int _Nscales)
{
	m = _m;
	Nscales = _Nscales;
	//	this->objMesh = objMesh;

	minPos = Vec3(0.0f, 0.0f, 0.0f);
	maxPos = Vec3(1.0f + Epsilon, 1.0f + Epsilon, 1.0f + Epsilon);
	cellSize = Vec3(0.1, 0.1, 0.1);

	ctLeaf = nullptr;
	ctGraph = nullptr;
	pcsOct = nullptr;
	kmeans = new KMeans;

#ifdef USE_SPARSE
	spLaplacian = nullptr;
	coefficients = new vector<T>;
#endif // USE_SPARSE

#ifdef SGWT_DEBUG
	fastSgwt = nullptr;
#endif //SGWT_DEBUG
}

PcsOctree::~PcsOctree()
{
	//不需要delete objMesh，objMesh指示借用别人创建好的对象。

	if (ctLeaf)
		delete ctLeaf;
	if (ctGraph)
		delete ctGraph;
	if (pcsOct)
		delete pcsOct;
	if (kmeans)
		delete kmeans;

#ifdef USE_SPARSE
	if (spLaplacian)
		delete	spLaplacian;

	if (coefficients)
		delete coefficients;

#ifdef SGWT_DEBUG
	if (fastSgwt)
		delete fastSgwt;
#endif //SGWT_DEBUG

#endif // USE_SPARSE
}

void PcsOctree::setParam(Vec3 min, Vec3 max, Vec3 cellSize)
{
	this->minPos = min;
	this->maxPos = max;
	this->cellSize = cellSize;
	if (!pcsOct)
	{
		pcsOct = new Octree<Node>(min, max, cellSize);
	}
	else
	{
		delete pcsOct;
		pcsOct = new Octree<Node>(min, max, cellSize);
	}
	if (ctLeaf)
		delete ctLeaf;
	if (ctGraph)
		delete ctGraph;
	if (pcsOct)
		delete pcsOct;

	ctLeaf = new CallTraverseGetInfoSetLeaf;
	pcsOct = new Octree<Node>(min, max, cellSize);

#ifdef USE_EIGEN
	ctGraph = new CallTGetGraph(&dMat, &weightAMat);
#endif //USE_EIGEN

#ifdef USE_SPARSE
	if (spLaplacian)
		delete	spLaplacian;
	ctGraph = new CallTGetGraph(coefficients);
#endif // USE_SPARSE

#ifdef USE_ARPACK
	ctGraph = new CallTGetGraph(&laplacianMat);
#endif //USE_ARPACK

	ctGraph->initParam(cellSize);
}

void PcsOctree::buildPcsOctFrmPC(ObjMesh * objeMesh)
{
	cout << endl << "Building octree...." << endl;
	// Initialize octree.
	// Minimum coordinates.
	// Add all particles to the octree.
	for (int v_it = 0; v_it < objeMesh->vertexList.size(); v_it++)
	{
		Node& n = pcsOct->getCell(objeMesh->vertexList[v_it]);

		//把点放到叶子节点所属的里面。
		n.pointPosList.push_back(objeMesh->vertexList[v_it]);
		n.colorList.push_back(objeMesh->colorList[v_it]);
		n.pointIdxList.push_back(v_it);
	}
	cout << "Building octree done." << endl;
	return;
}

void PcsOctree::getLeafboundary()
{
	cout << "Computing Leavf boundary...." << endl;

	// Prepare the callback class.
	pcsOct->traverse(ctLeaf);

#ifdef ZJW_DEBUG
	cout << "Computing Leavf boundary done!" << endl;
#endif // ZJW_DEUG
}

#ifdef RELATIVE_DIS_SIGNAL
void PcsOctree::changeNodePointsToRelativeDis()
{
	//遍历所有的叶子结点
	for (int leaf_it = 0; leaf_it < ctLeaf->nodeList.size(); leaf_it++)
	{
		//遍历这个叶子结点中的所有points
		Octree<Node>::OctreeNode * tempNode = ctLeaf->nodeList[leaf_it];
		
		for (int p_it = 0; p_it < tempNode->nodeData.pointPosList.size(); p_it++)
		{
			//减去这个Node边界上的最小值，就得到这个Point 在这个node中的相对位置。
#ifdef RELATIVE_DIS_SIGNAL
			//利用相对距离
			//减去这个Node边界上的最小值，就得到这个Point 在这个node中的相对位置。
			tempNode->nodeData.pointPosList[p_it] = tempNode->nodeData.pointPosList[p_it];
#elif
			//绝对距离
			//减去这个Node边界上的最小值，就得到这个Point 在这个node中的相对位置。
			tempNode->nodeData.pointPosList[p_it] = tempNode->nodeData.pointPosList[p_it] - tempNode->min;
#endif //relative dis signal
		}
	}
}
#endif //relative dis signal

void PcsOctree::initMat()
{
#ifdef ZJW_DEBUG
	cout << "leaf node num :  " << ctLeaf->nodeList.size() << endl;
	cout << "init the matriex !" << endl;
#endif
	nodeNum = ctLeaf->nodeList.size();
	if (nodeNum == 0)
	{
		cout << "CallTGetGraph::initParam error" << endl;
		return;
	}
#ifdef USE_EIGEN
	dMat = MatrixXd::Zero(nodeNum, nodeNum);
	weightAMat = MatrixXd::Zero(nodeNum, nodeNum);
	LaplacianMat = MatrixXd::Zero(nodeNum, nodeNum);
	eigenVecMat = MatrixXd::Zero(nodeNum, nodeNum);
	eigenValMat = MatrixXd::Zero(nodeNum, nodeNum);
#endif //USE_EIGEN

#ifdef USE_SPARSE
	spLaplacian = new SpMat(nodeNum, nodeNum);
#endif // USE_SPARSE

#ifdef USE_ARPACK
	laplacianMat.initParam(nodeNum);
#endif // USE_ARPACK
}

void PcsOctree::getGraphMat()
{
	initMat();

#ifdef ZJW_TIMER
	ZjwTimer timer;
	timer.Start();
	cout << "start get the D & weight matriex...." << endl;
#endif
	//设置实际的cell size 到ctGraph中
	ctGraph->initParam(ctLeaf->nodeList[0]->max - ctLeaf->nodeList[0]->min);
	cout << "real cell size :" << ctGraph->octCellSize.x << "  " <<ctGraph->octCellSize.y << "  " <<ctGraph->octCellSize.z<< endl;

	//遍历所有的叶子节点，找到每个叶子节点的相连的节点,并得到权重
	for (int i = 0; i < ctLeaf->nodeList.size(); i++)
	{
		//test
		/*if (i == 233)
		{
			int iasdf = 100;
		}*/
		//end test
		ctGraph->leafIdx = ctLeaf->nodeList[i]->leafFlag;
		//保存叶子节点中中心节点
		ctGraph->leafMidPoint = (ctLeaf->nodeList[i]->max + ctLeaf->nodeList[i]->min) / 2;
		pcsOct->traverse(ctGraph);
	}

#ifdef USE_EIGEN
	LaplacianMat = dMat - weightAMat;
#endif //use eigen

#ifdef  USE_SPARSE


	//生成稀疏的LaplacianMat
	assert(coefficients->size() > 0);
	spLaplacian->setFromTriplets(coefficients->begin(), coefficients->end());


	//test
	//timer.Start();
	//cout << "start get Mat Eigen Vertor and Value...." << endl;

	//EigenSolver<MatrixXd> es(*spLaplacian);
	//MatrixXd  eigenVecMat;
	//MatrixXd  eigenValMat;
	////way1
	////Matrix whose columns are the (possibly complex) eigenvectors.
	////es.eigenvectors().col(0)

	////way2
	////Const reference to matrix whose columns are the pseudo-eigenvectors.
	//eigenValMat = es.pseudoEigenvalueMatrix();
	//eigenVecMat = es.pseudoEigenvectors();


	//cout << "88888888888888888888888888888888888888888888888888888888888" << endl;
	//cout << "L :" << endl << *spLaplacian << endl;
	//cout << "The pseudo-eigenvalue matrix D is:" << endl << eigenValMat << endl;
	//cout << "The pseudo-eigenvector matrix V is:" << endl << eigenVecMat << endl;
	////cout << "Finally, V * D * V^(-1) = " << endl << eigenVecMat * eigenValMat * eigenVecMat.inverse() << endl;
	//cout << "88888888888888888888888888888888888888888888888888888888888" << endl;
	//cout << "getMatEigenVerValue time: " << timer.GetInMs() << " ms " << endl;
	//cout << "end get Mat Eigen Vertor and Value !" << endl;
	//timer.Stop();
	////end

#ifdef ZJW_PRINT_INFO
  //打印的是非0元素
	cout << "****************************" << endl;
	cout << "sparse mat lap: " << endl;
	for (int k = 0; k < spLaplacian->outerSize(); ++k)
	{
		for (SpMat::InnerIterator it(*spLaplacian, k); it; ++it)
		{
			//cout << it.value() << " " << it.row() << " " << it.col() << " " << it.index() << endl;   // row index
			cout << it.value() << " ";   // row index
		}
		cout << endl;
	}
#endif //ZJW_PRINT_INFO
#endif //  USE_SPARSE

#ifdef ZJW_TIMER
	timer.Stop();
	cout << "getGraphMat time: " << timer.GetInMs() << " ms " << endl;
	cout << "finish get the D & weight matriex !" << endl;
	//printMat();
#endif
}

void PcsOctree::getMatEigenVerValue()
{
#ifdef USE_EIGEN
	ZjwTimer timer;
	timer.Start();
	cout << "start get Mat Eigen Vertor and Value...." << endl;

	EigenSolver<MatrixXd> es(LaplacianMat);

	//way1
	//Matrix whose columns are the (possibly complex) eigenvectors.
	//es.eigenvectors().col(0)

	//way2
	//Const reference to matrix whose columns are the pseudo-eigenvectors.
	eigenValMat = es.pseudoEigenvalueMatrix();
	eigenVecMat = es.pseudoEigenvectors();

#ifdef ZJW_PRINT_INFO

	/*cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
	cout << "The matrix of eigenvectors, V, is:" << endl << es.eigenvectors() << endl << endl;*/

	cout << "88888888888888888888888888888888888888888888888888888888888" << endl;
	cout << "L :" << endl << LaplacianMat << endl;
	cout << "The pseudo-eigenvalue matrix D is:" << endl << eigenValMat << endl;
	cout << "The pseudo-eigenvector matrix V is:" << endl << eigenVecMat << endl;
	cout << "Finally, V * D * V^(-1) = " << endl << eigenVecMat * eigenValMat * eigenVecMat.inverse() << endl;
	cout << "88888888888888888888888888888888888888888888888888888888888" << endl;

#endif //ZJW_PRINT_INFO
	cout << "getMatEigenVerValue time: " << timer.GetInMs() << " ms " << endl;
	cout << "end get Mat Eigen Vertor and Value !" << endl;
	timer.Stop();

#endif // USE_EIGEN

#ifdef USE_ARPACK
	ZjwTimer timer;
	timer.Start();
	cout << "start get Mat Eigen Vertor and Value...." << endl;

	laplacianMat.dsaupdEvalsEvecs();

#ifdef  ZJW_PRINT_INFO
	laplacianMat.printMat();
	laplacianMat.printValueVector();
#endif //  ZJW_PRINT_INFO

	cout << "getMatEigenVerValue time: " << timer.GetInMs() << " ms " << endl;
	cout << "end get Mat Eigen Vertor and Value !" << endl;
	timer.Stop();

#endif //use_arpack
}

//拿到nodeidx这个点的two hop。拿到的点不包括自己，因为是无向图，所以这里的two hop不包括一步的
void PcsOctree::getTwoHopNeighborhood(int nodeIdx, set<int>* nodeList_out, SpMat * spLaplacian)
{
	assert(nodeIdx > -1 && nodeIdx < ctLeaf->nodeList.size());
	//只支持列主序
	assert(!spLaplacian->IsRowMajor);
	assert(nodeList_out);
	
	//保存第一跳的结果
	vector<int> tempList;
	//遍历这一列的所有非0元素
	for (SparseMatrix<double>::InnerIterator it(*spLaplacian, nodeIdx); it; ++it)
	{
		if (it.value() != 0)
		{
			tempList.push_back(it.row());
		}
	}

	//set<int> nodeList_out;
	nodeList_out->clear();
	for (int node_it = 0; node_it < tempList.size(); node_it++)
	{
		int index = tempList[node_it];
		for (SparseMatrix<double>::InnerIterator it(*spLaplacian, index); it; ++it)
		{
			if (it.value() != 0)
			{
				nodeList_out->insert(it.row());
			}
		}
	}
	//移除自己
	nodeList_out->erase(nodeIdx);
}

void PcsOctree::getTwoHopNeighborhoodWithOneStep(int nodeIdx, set<int>* nodeList_out, SpMat * spLaplacian)
{
	assert(!spLaplacian->IsRowMajor);
	assert(nodeList_out);

	//保存第一跳的结果
	vector<int> tempList;
	//遍历这一列的所有非0元素
	for (SparseMatrix<double>::InnerIterator it(*spLaplacian, nodeIdx); it; ++it)
	{
		if (it.value() != 0)
		{
			tempList.push_back(it.row());
		}
	}

	//set<int> nodeList_out;
	nodeList_out->clear();
	for (int node_it = 0; node_it < tempList.size(); node_it++)
	{
		nodeList_out->insert(tempList[node_it]);

		int index = tempList[node_it];
		for (SparseMatrix<double>::InnerIterator it(*spLaplacian, index); it; ++it)
		{
			if (it.value() != 0)
			{
				nodeList_out->insert(it.row());
			}
		}
	}
	//移除自己
	nodeList_out->erase(nodeIdx);
}

void PcsOctree::clearOct()
{
	if (pcsOct != nullptr)
		pcsOct->clear();
}

int PcsOctree::judege8Aeros(Vec3 & min, Vec3 & max, Vec3 & point)
{
	unsigned int  i = 0;
	Vec3 mid = 0.5 * (max + min);

	if (point.x >= mid.x)
	{
		i = i | 1;
	}

	if (point.y >= mid.y)
	{
		i = i | 2;
	}

	if (point.z >= mid.z)
	{
		i = i | 4;
	}

	return i;
}

int PcsOctree::judege8Aeros(Vec3 & mid, Vec3 & point)
{
	unsigned int  i = 0;

	if (point.x >= mid.x)
	{
		i = i | 1;
	}

	if (point.y >= mid.y)
	{
		i = i | 2;
	}

	if (point.z >= mid.z)
	{
		i = i | 4;
	}

	return i;
}

void PcsOctree::setPointTo8Areas()
{
#ifdef ZJW_DEBUG
	cout << "start set leafnode points to 8 areas ...." << endl;
#endif

	for (int leaf_it = 0; leaf_it < ctLeaf->nodeList.size(); leaf_it++)
	{
		Octree<Node>::OctreeNode* octNode = ctLeaf->nodeList[leaf_it];

		octNode->nodeData.leafNodePos8Areas.clear();
		octNode->nodeData.leafNodeColor8Areas.clear();

		octNode->nodeData.leafNodePos8Areas.resize(8);
		octNode->nodeData.leafNodeColor8Areas.resize(8);

		Vec3 mid = 0.5 * (octNode->max + octNode->min);
		//遍历该叶子节点上所有的点
		for (int p_it = 0; p_it < octNode->nodeData.pointPosList.size(); p_it++)
		{
#ifdef RELATIVE_DIS_SIGNAL
			//判断出该节点所属的象限
			int index = judege8Aeros(mid, octNode->nodeData.pointPosList[p_it] + octNode->min);

#else
			//绝对距离
			//判断出该节点所属的象限
			int index = judege8Aeros(mid, octNode->nodeData.pointPosList[p_it]);
#endif //relative dis signal
			//绝对距离的时候，放到是聚堆距离的信号。相对位置的时候，放的相对位置的信号。
			octNode->nodeData.leafNodePos8Areas[index].push_back((Vec3 *)(&octNode->nodeData.pointPosList[p_it]));
			octNode->nodeData.leafNodeColor8Areas[index].push_back((Vec3 *)(&octNode->nodeData.colorList[p_it]));
		}
	}

#ifdef ZJW_DEBUG
	cout << "end set leafnode points to 8 areas!!!!" << endl;
#endif
}

vector<VectorXd> PcsOctree::getSignalF(SignalType sType)
{
	//初始化信息
	vector<VectorXd> posSignal;
	VectorXd initF(nodeNum);
	initF.setZero();
	for (int i = 0; i < 8; i++)
	{
		posSignal.push_back(initF);
	}

	switch (sType)
	{
	case SignalX:
		//遍历所有的叶子节点
		for (int i = 0; i < ctLeaf->nodeList.size(); i++)
		{
			//遍历八个象限, 用对应象限的信号的值
			for (int j = 0; j < 8; j++)
			{
				//排除象限中没有点的情况
				if (ctLeaf->nodeList[i]->nodeData.pos8Flag[j])
				{
					posSignal[j](i) = ctLeaf->nodeList[i]->nodeData.pos8AreasSignal[j].x;
				}
			}
		}
		break;
	case SignalY:
		//遍历所有的叶子节点
		for (int i = 0; i < ctLeaf->nodeList.size(); i++)
		{
			//遍历八个象限
			for (int j = 0; j < 8; j++)
			{
				//排除象限中没有点的情况
				if (ctLeaf->nodeList[i]->nodeData.pos8Flag[j])
				{
					posSignal[j](i) = ctLeaf->nodeList[i]->nodeData.pos8AreasSignal[j].y;
				}
			}
		}
		break;
	case SignalZ:
		//遍历所有的叶子节点
		for (int i = 0; i < ctLeaf->nodeList.size(); i++)
		{
			//遍历八个象限
			for (int j = 0; j < 8; j++)
			{
				//排除象限中没有点的情况
				if (ctLeaf->nodeList[i]->nodeData.pos8Flag[j])
				{
					posSignal[j](i) = ctLeaf->nodeList[i]->nodeData.pos8AreasSignal[j].z;
				}
			}
		}
		break;
	case SignalR:
		//遍历所有的叶子节点
		for (int i = 0; i < ctLeaf->nodeList.size(); i++)
		{
			//遍历八个象限
			for (int j = 0; j < 8; j++)
			{
				//排除象限中没有点的情况
				if (ctLeaf->nodeList[i]->nodeData.pos8Flag[j])
				{
					posSignal[j](i) = ctLeaf->nodeList[i]->nodeData.color8AreasSignal[j].x;
				}
			}
		}
		break;
	case SignalG:
		//遍历所有的叶子节点
		for (int i = 0; i < ctLeaf->nodeList.size(); i++)
		{
			//遍历八个象限
			for (int j = 0; j < 8; j++)
			{
				//排除象限中没有点的情况
				if (ctLeaf->nodeList[i]->nodeData.pos8Flag[j])
				{
					posSignal[j](i) = ctLeaf->nodeList[i]->nodeData.color8AreasSignal[j].y;
				}
			}
		}
		break;
	case SignalB:
		//遍历所有的叶子节点
		for (int i = 0; i < ctLeaf->nodeList.size(); i++)
		{
			//遍历八个象限
			for (int j = 0; j < 8; j++)
			{
				//排除象限中没有点的情况
				if (ctLeaf->nodeList[i]->nodeData.pos8Flag[j])
				{
					posSignal[j](i) = ctLeaf->nodeList[i]->nodeData.color8AreasSignal[j].z;
				}
			}
		}
		break;
	default:
		break;
	}

	return posSignal;
}

bool PcsOctree::getAllSignalAndSaveSGWTCoeff()
{

#ifdef SGWT_DEBUG
	if (!fastSgwt)
	{
		assert(m > 0 && Nscales > 0);
		fastSgwt = new SgwtCheby(m, Nscales, *spLaplacian);
	}

#ifdef ZJW_TIMER
	ZjwTimer test;
	test.Start();
#endif //ZJW_TIMER

	//拿到这个信号，在所有象限中的信号，所有结点的信号。
	vector<SignalType> typeList;
	typeList.push_back(SignalX);
	typeList.push_back(SignalY);
	typeList.push_back(SignalZ);
	typeList.push_back(SignalR);
	typeList.push_back(SignalG);
	typeList.push_back(SignalB);

	int totalQuadrant = 8;

	//8个vectorXd， 分别表示不同象限的信号
	vector<VectorXd> fSignal;
	//遍历所有的信号
	for (int type_it = 0; type_it < typeList.size(); type_it++)
	{
		//8个vectorXd， 分别表示不同象限的信号。每个vector表示的是：所有在节点在这个象限中的该型号的具体的值
		fSignal = getSignalF(typeList[type_it]);

		//遍历所有的象限
		for (int quadrant_it = 0; quadrant_it < totalQuadrant; quadrant_it++)
		{
			//vector中表示的5个不同尺度下面所有顶点的向量
			(*fastSgwt).sgwt->saveSgwtCoeff(typeList[type_it], quadrant_it, &fSignal[quadrant_it], &((*fastSgwt).sgwt->coeff));
		}
	}

#ifdef ZJW_TIMER
	test.Stop();
	test.printTimeInMs("getAllSignalAndSaveSGWTCoeff time : ");
#endif //zjw_timer

#endif //SGWT_DEBUG
	return true;
}

//返回这个节点在特定信号下面所有象限信号，返回值应该是8个系数的向量，对应每个象限
bool PcsOctree::getSignalF(SignalType sType, int nodeIdx, VectorXd * sigVQ_out)
{
	//初始化信息
	if (sigVQ_out)
		delete sigVQ_out;
	sigVQ_out = new VectorXd(8);
	sigVQ_out->setZero();

	switch (sType)
	{
	case SignalX:
		//遍历八个象限, 用对应象限的信号的值
		for (int j = 0; j < 8; j++)
		{
			//排除象限中没有点的情况
			if (ctLeaf->nodeList[nodeIdx]->nodeData.pos8Flag[j])
			{
				(*sigVQ_out)(j) = ctLeaf->nodeList[nodeIdx]->nodeData.pos8AreasSignal[j].x;
			}
		}

		break;
	case SignalY:
		//遍历八个象限, 用对应象限的信号的值
		for (int j = 0; j < 8; j++)
		{
			//排除象限中没有点的情况
			if (ctLeaf->nodeList[nodeIdx]->nodeData.pos8Flag[j])
			{
				(*sigVQ_out)(j) = ctLeaf->nodeList[nodeIdx]->nodeData.pos8AreasSignal[j].y;
			}
		}
		break;
	case SignalZ:
		//遍历八个象限, 用对应象限的信号的值
		for (int j = 0; j < 8; j++)
		{
			//排除象限中没有点的情况
			if (ctLeaf->nodeList[nodeIdx]->nodeData.pos8Flag[j])
			{
				(*sigVQ_out)(j) = ctLeaf->nodeList[nodeIdx]->nodeData.pos8AreasSignal[j].z;
			}
		}
		break;
	case SignalR:
		//遍历八个象限, 用对应象限的信号的值
		for (int j = 0; j < 8; j++)
		{
			//排除象限中没有点的情况
			if (ctLeaf->nodeList[nodeIdx]->nodeData.pos8Flag[j])
			{
				(*sigVQ_out)(j) = ctLeaf->nodeList[nodeIdx]->nodeData.color8AreasSignal[j].x;
			}
		}
		break;
	case SignalG:
		for (int j = 0; j < 8; j++)
		{
			//排除象限中没有点的情况
			if (ctLeaf->nodeList[nodeIdx]->nodeData.pos8Flag[j])
			{
				(*sigVQ_out)(j) = ctLeaf->nodeList[nodeIdx]->nodeData.color8AreasSignal[j].y;
			}
		}
		break;
	case SignalB:
		for (int j = 0; j < 8; j++)
		{
			//排除象限中没有点的情况
			if (ctLeaf->nodeList[nodeIdx]->nodeData.pos8Flag[j])
			{
				(*sigVQ_out)(j) = ctLeaf->nodeList[nodeIdx]->nodeData.color8AreasSignal[j].z;
			}
		}
		break;
	default:
		break;
	}

	return true;
}

vector<VectorXd> PcsOctree::getSgwtCoeffWS(SignalType type, int quadrant)
{
	//拿到这个信号，在所有象限中的信号，所有结点的信号。
	vector<VectorXd> fSignal = getSignalF(type);

#ifdef SGWT_DEBUG
	if (!fastSgwt)
	{
		assert(m > 0 && Nscales > 0);
		fastSgwt = new SgwtCheby(m, Nscales, *spLaplacian);
	}
	//得到指定象限的，这个信号的 sgwt的系数

	/*ZjwTimer test;
	test.Start();*/
	vector<VectorXd> wf_s = (*fastSgwt)(fSignal[quadrant]);
	/*test.Stop();
	test.printTimeInMs("sgwt_cheby_op: get sgw coeff for signal f ");*/

#ifdef ZJW_DEBUG
	cout << "********* getSgwtCoeffWS *********" << endl;
	cout << "sgwt for signal " << type << " in quadrant " << quadrant << endl;
	fastSgwt->sgwt->getVectorVectorXdInfo(wf_s);
	cout << "********************************" << endl;
#endif //ZJW_DEUG

#endif //SGWT_DEBUG
	return wf_s;
}

bool PcsOctree::getFeatureVector(int nodeIdx, VectorXd *featureVector)
{
#ifdef SGWT_DEBUG
	if (!fastSgwt)
	{
		assert(m > 0 && Nscales > 0);
		fastSgwt = new SgwtCheby(m, Nscales, *spLaplacian);
	}

#ifdef ZJW_TIMER
	/*ZjwTimer test;
	test.Start();*/
#endif //ZJW_TIMER

	//拿到这个信号，在所有象限中的信号，所有结点的信号。
	vector<SignalType> typeList;
	typeList.push_back(SignalX);
	typeList.push_back(SignalY);
	typeList.push_back(SignalZ);
	typeList.push_back(SignalR);
	typeList.push_back(SignalG);
	typeList.push_back(SignalB);

	//记录g func的scale的个数，在加上一个h func
	int totalScale = (fastSgwt->sgwt->t.size() + 1);
	int totalSignal = typeList.size();
	int totalQuadrant = 8;
	featureVector->resize(totalQuadrant * totalSignal * totalScale);

	//向量下标
	int idx = -1;
	//8个vectorXd， 分别表示不同象限的信号
	vector<VectorXd> fSignal;
	//遍历所有的信号
	for (int type_it = 0; type_it < typeList.size(); type_it++)
	{
		//8个vectorXd， 分别表示不同象限的信号。每个vector表示的是：所有在节点在这个象限中的该型号的具体的值
		fSignal = getSignalF(typeList[type_it]);

		//遍历所有的象限
		for (int quadrant_it = 0; quadrant_it < 8; quadrant_it++)
		{
			//vecotr中表示的5个不同尺度下面所有顶点的向量
			vector<VectorXd> wf_s = (*fastSgwt)(fSignal[quadrant_it]);

			//得到该节点在五个尺度下的系数：
			for (int s_it = 0; s_it < totalScale; s_it++)
			{
				idx = type_it * totalQuadrant * totalScale + quadrant_it * totalScale + s_it;
				(*featureVector)(idx) = wf_s[s_it](nodeIdx);
			}
		}
	}

#ifdef ZJW_TIMER
	/*test.Stop();
	cout << nodeIdx;
	test.printTimeInMs(" get feature vector time : ");*/
#endif //zjw_timer

#endif //SGWT_DEBUG
	return true;
}

bool PcsOctree::getFeatureVector2(int nodeIdx, VectorXd * featureVector)
{
#ifdef SGWT_DEBUG
	if (!fastSgwt)
	{
		assert(m > 0 && Nscales > 0);
		fastSgwt = new SgwtCheby(m, Nscales, *spLaplacian);
	}

	//#ifdef ZJW_TIMER
	//	ZjwTimer test;
	//	test.Start();
	//#endif //ZJW_TIMER

		//拿到这个信号，在所有象限中的信号，所有结点的信号。
	vector<SignalType> typeList;
	typeList.push_back(SignalX);
	typeList.push_back(SignalY);
	typeList.push_back(SignalZ);
	typeList.push_back(SignalR);
	typeList.push_back(SignalG);
	typeList.push_back(SignalB);

	int totalScale = (fastSgwt->sgwt->t.size() + 1);
	int totalSignal = typeList.size();
	int totalQuadrant = 8;
	featureVector->resize(totalQuadrant * totalSignal * totalScale);

	//向量下标
	int idx = -1;
	//表示这个顶点在这个信号，这个象限的系数 维度：5*1
	VectorXd w_s;
	//遍历所有的信号
	for (int type_it = 0; type_it < typeList.size(); type_it++)
	{
		//遍历所有的象限
		for (int quadrant_it = 0; quadrant_it < 8; quadrant_it++)
		{
			//wf_s保存的是： vecotr中表示的5个不同尺度下面所有顶点的的系数
			(*fastSgwt)(nodeIdx, typeList[type_it], quadrant_it, &w_s);
			//得到该节点在五个尺度下的系数：
			for (int s_it = 0; s_it < totalScale; s_it++)
			{
				//int nodeIdx, int signalType, int quadrantType, VectorXd* sgwt_out
				idx = type_it * totalQuadrant * totalScale + quadrant_it * totalScale + s_it;
				(*featureVector)(idx) = w_s(s_it);
			}
		}
	}

	//#ifdef ZJW_TIMER
	//	test.Stop();
	//	cout << nodeIdx;
	//	test.printTimeInMs(" getFeatureVector2  time : ");
	//#endif //zjw_timer

#endif //SGWT_DEBUG

	return true;
}

void PcsOctree::getSgwtCoeffWS()
{
	posSignalX = getSignalF(SignalType::SignalX);
	posSignalY = getSignalF(SignalType::SignalY);
	posSignalZ = getSignalF(SignalType::SignalZ);

#ifdef SGWT_DEBUG
	if (!fastSgwt)
	{
		assert(m > 0 && Nscales > 0);
		fastSgwt = new SgwtCheby(m, Nscales, *spLaplacian);
	}
#endif //SGWT_DEBUG

	//得到信号X 在8个象限中的信号
	for (int i = 0; i < 8; i++)
	{
		//单个象限，单个信号下的，wf
#ifdef SGWT_DEBUG
		vector<VectorXd> wf_s = (*fastSgwt)(posSignalX[i]);
#ifdef ZJW_DEBUG
		fastSgwt->sgwt->getVectorVectorXdInfo(wf_s);
		//fastSgwt->sgwt->printVectorVectorXd(wf_s);
#endif //ZJW_DEUG
#endif //SGWT_DEBUG
	}
}

void PcsOctree::doKmeans(int clusterNum)
{
#ifdef ZJW_TIMER
	cout << "start do Kmeans ...." << endl;
	cout << "kmeans cluster set to " << clusterNum << endl;
	ZjwTimer timer2;
	timer2.Start();
#endif

	//设置分类的组数
	kmeans->setClusterNum(clusterNum);
	//k means进行聚类
	kmeans->cluster(ctLeaf->midVList);

	//结果保存在clusterRes，保存了叶子节点的序号
	//dosomething

#ifdef ZJW_DEBUG
	//打印分类信息
	cout << "kmeans result: " << endl;
	for (int c_it = 0; c_it < kmeans->clusterRes.size(); c_it++)
	{
		cout << "cluster " << c_it << " : " << endl;

		for (int n_it = 0; n_it < kmeans->clusterRes[c_it].size(); n_it++)
		{
			cout << kmeans->clusterRes[c_it][n_it]<<" ";
		}
		cout <<endl;
	}
#endif // ZJW_DEBUG

#ifdef ZJW_TIMER
	timer2.Stop();
	timer2.printTimeInMs("end do Kmeans: ");
#endif
}

int PcsOctree::judegePointToLeafNode(Vec3 * point, int & idx)
{
	AddPoints judgeP(*point);
	pcsOct->traverse(&judgeP);
	idx = judgeP.nodeIdx;

	//如果没有匹配成功，返回的是-1.
	return idx;
}

void PcsOctree::getLeafSignal()
{
#ifdef ZJW_DEBUG
	cout << "start get Leaf Signal ...." << endl;
#endif

	//遍历每个叶子节点
	for (int leaf_it = 0; leaf_it < ctLeaf->nodeList.size(); leaf_it++)
	{
		Octree<Node>::OctreeNode* octNode = ctLeaf->nodeList[leaf_it];

		//遍历这个node上所有的点，保存这个node的平均颜色
		Color nodeTotalColor(0,0,0);
		for (int i = 0; i < octNode->nodeData.colorList.size(); i++)
		{
			nodeTotalColor += octNode->nodeData.colorList[i];
		}
		ctLeaf->nodeColorList.push_back(nodeTotalColor / octNode->nodeData.colorList.size());

		//遍历当前叶子节点的八个象限,分别得到每个象限的颜色值以及位置信息
		for (int i = 0; i < octNode->nodeData.leafNodePos8Areas.size(); i++)
		{
			//这个象限没有point
			if (octNode->nodeData.leafNodePos8Areas[i].size() == 0)
			{
				octNode->nodeData.pos8Flag.push_back(false);
				octNode->nodeData.pos8AreasSignal.push_back(Vec3(0, 0, 0));
				octNode->nodeData.color8AreasSignal.push_back(Vec3(0, 0, 0));
				continue;
			}

			octNode->nodeData.pos8Flag.push_back(true);

			//遍历这个象限的所有point,来计算信号。下面用的是求平均的方法，来计算每个象限的所有点的平均信号
			Vec3 posSignal(0, 0, 0);
			Vec3 colorSignal(0, 0, 0);

			for (int p_it = 0; p_it < octNode->nodeData.leafNodePos8Areas[i].size(); p_it++)
			{
				posSignal += *(octNode->nodeData.leafNodePos8Areas[i][p_it]);
				colorSignal += *(octNode->nodeData.leafNodeColor8Areas[i][p_it]);
			}
			posSignal /= octNode->nodeData.leafNodePos8Areas[i].size();
			colorSignal /= octNode->nodeData.leafNodeColor8Areas[i].size();

			octNode->nodeData.pos8AreasSignal.push_back(posSignal);
			octNode->nodeData.color8AreasSignal.push_back(colorSignal);
		}
	}

#ifdef ZJW_DEBUG
	cout << "end get Leaf Signal!!!!" << endl;
#endif
}

vector<Vec3>& PcsOctree::getLeafAvgColor() const
{
	// TODO: 在此处插入 return 语句
	return ctLeaf->nodeColorList;
}

void PcsOctree::printMat()
{
	cout << "****************************" << endl;
	cout << "node num : " << nodeNum << endl;
	cout << "****************************" << endl;

#ifdef USE_EIGEN
	cout << "D mat : " << endl;
	cout << dMat << endl;
	cout << "****************************" << endl;
	cout << "weight mat: " << endl;
	cout << weightAMat << endl;
#endif // USE_EIGEN
}

//------------------------------------------------------------------

CallTraverseGetInfoSetLeaf::CallTraverseGetInfoSetLeaf()
{
	minVList.clear();
	maxVList.clear();
	midVList = new vector<Vec3>;
	leafIncr = 0;
}

CallTraverseGetInfoSetLeaf::~CallTraverseGetInfoSetLeaf()
{
	if (midVList)
		delete midVList;

	for (int i = 0; i < midVList->size(); i++)
	{
		delete (*midVList)[i];
	}
}

bool CallTraverseGetInfoSetLeaf::operator()(const Vec3 min, const Vec3 max, Octree<Node>::OctreeNode * currNode)
{
	bool flag = true;
	//判断是会否是叶子节点
	if (currNode->children[0] || currNode->children[1] || currNode->children[2] || currNode->children[3]
		|| currNode->children[4] || currNode->children[5] || currNode->children[6] || currNode->children[7])
	{
		flag = true;
	}
	else
	{
		//是叶子节点，那么保留叶子节点的矩形范围，然后退出
		flag = false;
		minVList.push_back(currNode->min);
		maxVList.push_back(currNode->max);

		//保存叶子节点的中间
		Vec3 temp = (currNode->min + currNode->max) / 2;
		midVList->push_back(temp);
		nodeList.push_back(currNode);
		currNode->leafFlag = leafIncr;
		leafIncr++;
	}
	return flag;
}

//------------------------------------------------------------------

//CallTGetGraph::CallTGetGraph(Vec3 octreeCellSize)
//{
//	this->octCellSize = octreeCellSize;
//}

CallTGetGraph::~CallTGetGraph()
{
#ifdef USE_EIGEN
	delete dMatPtr;
	delete weightAMatPtr;
#endif // USE_EIGEN

#ifdef USE_SPARSE
	if (spLap)
		spLap = nullptr;
	//因为coeff指向的是其他有用的地址，不能直接delete,而是赋值为空，因为delete会把指针指向的内存块也删除掉
	if (coeff)
		coeff = nullptr;

#endif // USE_SPARSE

#ifdef  USE_ARPACK
	delete laplacianMat;
#endif //  use_arpack
}
void CallTGetGraph::initParam(Vec3 & octreeCellSize)
{
	this->octCellSize = octreeCellSize;
}
//返回false,终止递归遍历;返回true,继续递归遍历子节点
bool CallTGetGraph::operator()(const Vec3 min, const Vec3 max, Octree<Node>::OctreeNode * currNode)
{
	double curNodeLegth = Length(max - min);
	//根据碰撞检测，来考虑这个值
	double maxDist = curNodeLegth / 2 + Length(octCellSize) / 2 + Epsilon;

	//-----------快速拒绝---------
	Vec3 midPoint = (currNode->min + currNode->max) / 2;
	Vec3 delta = midPoint - leafMidPoint;
	//太远了，不需要遍历子节点
	if (Length(delta) > maxDist)
		return false;

	//-------------如果可能找到邻居--------------

	//不是叶子节点，需要遍历子节点，继续进行快速拒绝
	int idx = currNode->leafFlag;
	if (idx == -1)
		return true;

	//当前节点是叶子节点（因为是无向图，是对称矩阵，只需要考虑一边）
	if (idx >= leafIdx)
		return false;

	//得到矩阵
#ifdef  USE_EIGEN
	(*weightAMatPtr)(idx, leafIdx) = 1 / Length(delta);
	(*weightAMatPtr)(leafIdx, idx) = 1 / Length(delta);
	(*dMatPtr)(leafIdx, leafIdx) += 1 / Length(delta);
	(*dMatPtr)(idx, idx) += 1 / Length(delta);
#endif //  USE_EIGEN

#ifdef USE_SPARSE
	//L  = D-W 下面直接保存为L了
	//test
	//if (idx == 233 || leafIdx == 233)
	//{
	//	if (idx == 233)
	//		cout << "node " << leafIdx;
	//	else
	//		cout << "node " << idx;
	//	//cout << "cellsize(斜对角 0.1*sgrt(3))： " << Length(octCellSize) << endl;
	//	cout << "相邻结点中点之间的距离：" << Length(delta) << endl;
	//}
	//end test
	coeff->push_back(T(idx, leafIdx, -1 / Length(delta)));
	coeff->push_back(T(leafIdx, idx, -1 / Length(delta)));
	coeff->push_back(T(leafIdx, leafIdx, 1 / Length(delta)));
	coeff->push_back(T(idx, idx, 1 / Length(delta)));
#endif // USE_SPARSE

#ifdef USE_ARPACK
	//L  = D-W 下面直接保存为L了

	//-W
	laplacianMat->setMat(idx, leafIdx, -1 / Length(delta));
	laplacianMat->setMat(leafIdx, idx, -1 / Length(delta));

	//D
	laplacianMat->setAddMat(leafIdx, 1 / Length(delta));
	laplacianMat->setAddMat(idx, 1 / Length(delta));
#endif // USE_ARPACK

	return true;
}

AddPoints::AddPoints(Vec3 point)
{
	nodeIdx = -1;
	this->point = point;
}

AddPoints::~AddPoints()
{
}


//返回false,终止递归遍历;返回true,继续递归遍历子节点
bool AddPoints::operator()(const Vec3 min, const Vec3 max, Octree<Node>::OctreeNode * currNode)
{
	//不在这个node里面，不需要遍历子节点
	if (!(point >= min && point < max))
		return false;

	//-------------如果在这个Node里面--------------

	//不是叶子节点，需要遍历子节点，继续进行快速拒绝
	int idx = currNode->leafFlag;
	if (idx == -1)
		return true;

	//当前节点是叶子节点，得到idx
	nodeIdx = idx;
	return true;
}