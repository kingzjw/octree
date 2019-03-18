#include "zjw_pcs_double_buffer_octree.h"

void PcsDBufferOctree::setTargetFrameId(int id)
{
	assert(id == 0 || id == 1);
	//一个是0，一个是1
	targetFrameId = id;
	refrenceFrameId = (id + 1) % 2;

	pcsDBufferOct->setTargetFrameId(id);
}

void PcsDBufferOctree::setReferenceFrameId(int id)
{
	assert(id == 0 || id == 1);
	//一个是0，一个是1
	targetFrameId = id;
	refrenceFrameId = (id + 1) % 2;

	pcsDBufferOct->setReferenceFrameId(id);
}

PcsDBufferOctree::PcsDBufferOctree(Vec3 minPos, Vec3 maxPosVec3, Vec3 cellSize)
{
	refre_leaf_count_=0;
	target_leaf_count_=0;
	refFramePointNum = 0;
	targetFramePointNum = 0;

	this->minPos = minPos;
	this->maxPos = maxPosVec3;
	this->cellSize = cellSize;

	targetFrameLoaded = false;
	referenceFrameLoaded = false;

	targetFrameId = 0;
	refrenceFrameId = 1;
	pcsDBufferOct = new DoubleBufferOctree<DBufferNodeData>(minPos, maxPosVec3, cellSize, targetFrameId, refrenceFrameId);
	ctLeaf = new TraverseGetInfoSetLeaf();
}

PcsDBufferOctree::PcsDBufferOctree(Vec3 minPos, Vec3 maxPosVec3, double cellSize)
{
	refre_leaf_count_ = 0;
	target_leaf_count_ = 0;
	
	ctLeaf = new TraverseGetInfoSetLeaf();


	refFramePointNum = 0;
	targetFramePointNum = 0;

	this->minPos = minPos;
	this->maxPos = maxPosVec3;
	this->cellSize.x = cellSize;
	this->cellSize.y = cellSize;
	this->cellSize.z = cellSize;

	targetFrameLoaded = false;
	referenceFrameLoaded = false;

	targetFrameId = 0;
	refrenceFrameId = 1;
	pcsDBufferOct = new DoubleBufferOctree<DBufferNodeData>(minPos, maxPosVec3, cellSize, targetFrameId, refrenceFrameId);
}

PcsDBufferOctree::PcsDBufferOctree()
{
	refre_leaf_count_ = 0;
	target_leaf_count_ = 0;

	ctLeaf = new TraverseGetInfoSetLeaf();
	refFramePointNum = 0;
	targetFramePointNum = 0;

	this->minPos = Vec3(0, 0, 0);
	this->maxPos = Vec3(1, 1, 1);
	this->cellSize.x = 0.5;
	this->cellSize.y = 0.5;
	this->cellSize.z = 0.5;

	targetFrameLoaded = false;
	referenceFrameLoaded = false;

	targetFrameId = 0;
	refrenceFrameId = 1;
	pcsDBufferOct = new DoubleBufferOctree<DBufferNodeData>(minPos, maxPos, cellSize, targetFrameId, refrenceFrameId);
}

PcsDBufferOctree::~PcsDBufferOctree()
{
	if (pcsDBufferOct)
		delete pcsDBufferOct;
	if (ctLeaf)
		delete ctLeaf;
}

void PcsDBufferOctree::buildDBufferOctree(ObjMesh * referenceMesh, ObjMesh * targetMesh)
{
	assert(!targetFrameLoaded && !referenceFrameLoaded);

	buildDBufferOctree(false, referenceMesh);
	buildDBufferOctree(true, targetMesh);
}

void PcsDBufferOctree::buildDBufferOctree(vector<Vertex> & referenceMeshVerList, vector<Vertex> & targetMeshVerList)
{
	assert(!targetFrameLoaded && !referenceFrameLoaded);

	buildDBufferOctree(false, referenceMeshVerList);
	buildDBufferOctree(true, targetMeshVerList);
}

uint32_t PcsDBufferOctree::getFramePointNum(bool isTarget)
{
	if (isTarget)
	{
		return targetFramePointNum;
	}
	else
	{
		return refFramePointNum;
	}
}

void PcsDBufferOctree::getBoundingBox(double & min_x, double & min_y, double & min_z, double & max_x, double & max_y, double & max_z)
{
	min_x = minPos.x;
	min_y = minPos.y;
	min_z = minPos.z;
	max_x = maxPos.x;
	max_y = maxPos.y;
	max_z = maxPos.z;
}

void PcsDBufferOctree::setBoundingBox(const double &min_x, const double &min_y, const double &min_z,
	const double & max_x, const double & max_y, const double &max_z)
{
	minPos.x = min_x;
	minPos.y = min_y;
	minPos.z = min_z;

	maxPos.x = max_x;
	maxPos.y = max_y;
	maxPos.z = max_z;

	pcsDBufferOct->min = this->minPos;
	pcsDBufferOct->max = this->maxPos;
}

void PcsDBufferOctree::setCellSize(const double cellSize)
{
	this->cellSize.x = cellSize;
	this->cellSize.y = cellSize;
	this->cellSize.z = cellSize;
	pcsDBufferOct->cellSize = this->cellSize;
}

void PcsDBufferOctree::getLeafboundaryAndLeafNode(bool isTarget)
{
	ctLeaf->initParam();
	ctLeaf->setIsTarget(isTarget);

	// Prepare the callback class.
	pcsDBufferOct->traverse(ctLeaf);

	if (isTarget)
	{
		target_leaf_count_ = ctLeaf->nodeList.size();
	}
	else
	{
		refre_leaf_count_ = ctLeaf->nodeList.size();
	}
}

void PcsDBufferOctree::buildDBufferOctree(bool isTarget, ObjMesh * objMesh)
{
	//如果已经加载过了，那么需要现清理
	if (isTarget && targetFrameLoaded)
	{
		clearOneFrame(isTarget);
	}
	else if ((!isTarget) && referenceFrameLoaded)
	{
		clearOneFrame(isTarget);
	}

	if (isTarget)
	{
		targetMeshVerList = &objMesh->vertexList;
		targetFramePointNum = objMesh->vertexList.size();
		targetFrameLoaded = true;
	}
	else
	{
		referenceMeshVerList = &objMesh->vertexList;
		refFramePointNum = objMesh->vertexList.size();
		referenceFrameLoaded = true;
	}

	//现在加载Frame
	for (int v_it = 0; v_it < objMesh->vertexList.size(); v_it++)
	{
		//利用深度优先，创建八叉树
		DBufferNodeData * node = pcsDBufferOct->getCell(isTarget, objMesh->vertexList[v_it]);

		//把点放到叶子节点所属的里面
		node->pointPosList.push_back(objMesh->vertexList[v_it]);
		node->pointIdxList.push_back(v_it);
	}
}

void PcsDBufferOctree::buildDBufferOctree(bool isTarget, vector<Vertex>& vertexList)
{
	//如果已经加载过了，那么需要现清理
	if (isTarget && targetFrameLoaded)
	{
		clearOneFrame(isTarget);
	}
	else if ((!isTarget) && referenceFrameLoaded)
	{
		clearOneFrame(isTarget);
	}

	if (isTarget)
	{
		targetMeshVerList = &vertexList;
		targetFramePointNum = vertexList.size();
		targetFrameLoaded = true;
	}
	else
	{
		referenceMeshVerList = &vertexList;
		refFramePointNum = vertexList.size();
		referenceFrameLoaded = true;
	}

	//现在加载Frame
	for (int v_it = 0; v_it < vertexList.size(); v_it++)
	{
		//利用深度优先，创建八叉树
		DBufferNodeData *node = pcsDBufferOct->getCell(isTarget, vertexList[v_it]);

		//把点放到叶子节点所属的里面。
		node->pointPosList.push_back(vertexList[v_it]);
		node->pointIdxList.push_back(v_it);
	}
}

void PcsDBufferOctree::clearOneFrame(bool isTarget)
{
	if (isTarget)
	{
		targetFrameLoaded = false;
		targetFramePointNum = 0;
		target_leaf_count_=0;
	}
	else
	{
		referenceFrameLoaded = false;
		refFramePointNum = 0;
		refre_leaf_count_=0;
	}

	//遍历所有的 结点，从下到上，删除无用的节点，或者是free NodeData(这里还要改写flag)
	TraverseAllNodeForFree freeTraverse(isTarget);
	pcsDBufferOct->traverseDepthFirst(&freeTraverse);
}

bool PcsDBufferOctree::getByteSteamOfOctreeXOR(vector<char>& byteList_out)
{
	byteList_out.clear();

	assert(targetFrameLoaded && referenceFrameLoaded);
	vector<OctreeDoubelBufferNode<DBufferNodeData>*>  allNodeList_out;
	pcsDBufferOct->traverseBreathFirst(allNodeList_out);

	for (int i = 0; i < allNodeList_out.size(); i++)
	{
		byteList_out.push_back(allNodeList_out[i]->getChildByteXOR());
	}
	return true;
}

bool PcsDBufferOctree::getByteSteamOfOctreeSigleFrame(vector<char>& byteList_out, bool isTarget)
{
	byteList_out.clear();
	vector<OctreeDoubelBufferNode<DBufferNodeData>*>  allNodeList_out;

	if (isTarget)
	{
		assert(targetFrameLoaded);
		pcsDBufferOct->traverseBreathFirstSigleFrame(allNodeList_out,true);
	}
	else
	{
		assert(referenceFrameLoaded);
		pcsDBufferOct->traverseBreathFirstSigleFrame(allNodeList_out, false);
	}

	for (int i = 0; i < allNodeList_out.size(); i++)
	{
		byteList_out.push_back(allNodeList_out[i]->getChildByteSigleFrame(isTarget));
	}

	return true;
}

void PcsDBufferOctree::clearDBufferOctree()
{
	//clear
	if (pcsDBufferOct)
		delete pcsDBufferOct;

	//init

	refFramePointNum = 0;
	targetFramePointNum = 0;
	targetFrameLoaded = false;
	referenceFrameLoaded = false;
	targetFrameId = 0;
	refrenceFrameId = 1;

	pcsDBufferOct = new DoubleBufferOctree<DBufferNodeData>(minPos, maxPos, cellSize, targetFrameId, refrenceFrameId);
}

bool PcsDBufferOctree::recoverySigleFrameInDBufferOctree(vector<char>& byteListForSigleFrame, bool isTarget)
{
	assert(!byteListForSigleFrame.empty());
	//init

	//以防万一，清空double buffer octree
	clearDBufferOctree();

	referenceFrameLoaded = true;

	//创建根节点上的关于isTarget这帧的Nodedata
	pcsDBufferOct->root->setNodeData(isTarget);

	//according to the byteList to recovery the target frame
	std::queue<OctreeDoubelBufferNode<DBufferNodeData>*> nodeQue;
	nodeQue.push(pcsDBufferOct->root);
	int c_it = 0;
	while (!nodeQue.empty())
	{
		OctreeDoubelBufferNode<DBufferNodeData>* tempNode = nodeQue.front();
		pcsDBufferOct->recoveryDBufferOctreeNodeForSigleFrame(nodeQue.front(), byteListForSigleFrame[c_it++], false);
		nodeQue.pop();
		for (int n_it = 0; n_it < 8; n_it++)
		{
			//孩子节点不为空的，放入到队列中
			if (tempNode->getChildren(n_it))
			{
				nodeQue.push(tempNode->getChildren(n_it));
			}
		}
	}
	return true;
}

bool PcsDBufferOctree::recoverySecondFrameInDBufferOctreeXOR(vector<char>& byteStreamXOR, bool isTarget)
{
	assert(!byteStreamXOR.empty());
	if (isTarget)
	{
		assert(referenceFrameLoaded && !targetFrameLoaded);
		targetFrameLoaded = true;

	}
	else
	{
		assert(targetFrameLoaded && !referenceFrameLoaded);
		referenceFrameLoaded = true;
	}

	//创建根节点上的关于isTarget这帧的Nodedata
	pcsDBufferOct->root->setNodeData(isTarget);

	//according to the byteList to recovery the target frame
	std::queue<OctreeDoubelBufferNode<DBufferNodeData>*> nodeQue;
	nodeQue.push(pcsDBufferOct->root);
	int c_it = 0;
	while (!nodeQue.empty())
	{
		OctreeDoubelBufferNode<DBufferNodeData>* tempNode = nodeQue.front();
		pcsDBufferOct->recoveryDBufferOctreeNodeXOR(nodeQue.front(), byteStreamXOR[c_it++]);
		nodeQue.pop();
		for (int n_it = 0; n_it < 8; n_it++)
		{
			//孩子节点不为空的，放入到队列中
			if (tempNode->getChildren(n_it))
			{
				nodeQue.push(tempNode->getChildren(n_it));
			}
		}
	}
	return false;
}

bool PcsDBufferOctree::recoveryDBufferOctreeForTarget(vector<Vertex>& swapFrameVerList, vector<char>& byteList)
{
	assert(!byteList.empty());
	//init
	clearDBufferOctree();
	targetFrameLoaded = true;

	//build swap frame octree in referenceframe in double buffer octree
	buildDBufferOctree(false, swapFrameVerList);

	//according to the byteList to recovery the target frame
	std::queue<OctreeDoubelBufferNode<DBufferNodeData>*> nodeQue;
	nodeQue.push(pcsDBufferOct->root);
	int c_it = 0;
	while (!nodeQue.empty())
	{
		OctreeDoubelBufferNode<DBufferNodeData>* tempNode = nodeQue.front();
		pcsDBufferOct->recoveryDBufferOctreeNodeXOR(nodeQue.front(), byteList[c_it++]);
		nodeQue.pop();
		for (int n_it = 0; n_it < 8; n_it++)
		{
			//孩子节点不为空的，放入到队列中
			if (tempNode->getChildren(n_it))
			{
				nodeQue.push(tempNode->getChildren(n_it));
			}
		}
	}
	return true;
}

bool PcsDBufferOctree::getTargetOctreeByClearReference()
{
	assert(targetFrameLoaded);

	if (referenceFrameLoaded)
		clearOneFrame(false);

	return true;
}


void PcsDBufferOctree::getNextTwoFrame(ObjMesh * newObjMesh)
{
	assert(targetFrameLoaded && referenceFrameLoaded);
	//交换target frame 和reference frame，因为target 马上要变成了reference，旧的reference没有用了。
	if (targetFrameId == 0)
		setTargetFrameId(1);
	else
		setTargetFrameId(0);

	//删除此时的target frame
	clearOneFrame(true);

	//加载新的target frame
	buildDBufferOctree(true, newObjMesh);
}

void PcsDBufferOctree::getNextTwoFrame(vector<Vertex>& newVertexList)
{
	assert(targetFrameLoaded && referenceFrameLoaded);
	//交换target frame 和reference frame，因为target 马上要变成了reference，旧的reference没有用了。
	if (targetFrameId == 0)
		setTargetFrameId(1);
	else
		setTargetFrameId(0);

	//删除此时的target frame
	clearOneFrame(true);

	//加载新的target frame
	buildDBufferOctree(true, newVertexList);

}

TraverseAllNodeForFree::TraverseAllNodeForFree(bool isTarget)
{
	this->isTarget = isTarget;
}

TraverseAllNodeForFree::~TraverseAllNodeForFree()
{

}

bool TraverseAllNodeForFree::operator()(const Vec3 min, const Vec3 max, OctreeDoubelBufferNode<DBufferNodeData>* currNode)
{
	//如果是叶子节点，那么终止。叶子节点由上一层来处理
	if (currNode->judegeLeafNode())
		return true;

	//清理这个node上的子节点的无用的数据
	if (isTarget)
	{
		//判断八个子节点
		for (int i = 0; i < 8; i++)
		{
			//这个子Node上 是否有target的数据
			if (currNode->flag[currNode->getTargetFlagIdx()][i])
			{
				//判断这个子节点上有没有reference上的数据
				if (currNode->flag[currNode->getReferenceFlagIdx()][i])
				{
					//有reference上的数据，那么只要重制flag信息，然后释放NodeData指针
					currNode->flag[currNode->getTargetFlagIdx()][i] = false;
					delete currNode->getTargetFrameNodeDataPtr();
				}
				else
				{
					//那么直接释放整个子节点
					delete  currNode->getChildren(i);
				}
			}
		}
	}
	else
	{
		//删除reference   判断八个子节点
		for (int i = 0; i < 8; i++)
		{
			//这个子Node上 是否有target的数据
			if (currNode->flag[currNode->getReferenceFlagIdx()][i])
			{
				//判断这个子节点上有没有reference上的数据
				if (currNode->flag[currNode->getTargetFlagIdx()][i])
				{
					//有reference上的数据，那么只要重制flag信息，然后释放NodeData指针
					currNode->flag[currNode->getReferenceFlagIdx()][i] = false;
					delete currNode->getRefrenceFrameNodeDataPtr();
				}
				else
				{
					//那么直接释放整个子节点
					delete  currNode->getChildren(i);
				}
			}
		}
	}
	return true;
}

TraverseGetInfoSetLeaf::TraverseGetInfoSetLeaf(bool isTarget_in)
{
	minVList.clear();
	maxVList.clear();
	nodeList.clear();
	leafIncr = 0;
	isTarget = isTarget_in;
}

TraverseGetInfoSetLeaf::~TraverseGetInfoSetLeaf()
{
	//nodelist中内容的指针是Node在管理，这里不用删除
}

void TraverseGetInfoSetLeaf::setIsTarget(bool flag)
{
	this->isTarget = flag;
}

void TraverseGetInfoSetLeaf::initParam()
{
	minVList.clear();
	maxVList.clear();
	nodeList.clear();
	leafIncr = 0;
}

bool TraverseGetInfoSetLeaf::operator()(const Vec3 min, const Vec3 max, OctreeDoubelBufferNode<DBufferNodeData>* currNode)
{
	bool flag = true;
	//判断是会否是叶子节点
	/*if (currNode->getChildren(0) || currNode->getChildren(1) || currNode->getChildren(2) || currNode->getChildren(3)
		|| currNode->getChildren(4) || currNode->getChildren(5) || currNode->getChildren(6) || currNode->getChildren(7))*/

	int idx = -1;
	if (isTarget)
		idx = currNode->getTargetFlagIdx();
	else
		idx = currNode->getReferenceFlagIdx();

	if (currNode->flag[idx][0] || currNode->flag[idx][1] || currNode->flag[idx][2] || currNode->flag[idx][3]
		|| currNode->flag[idx][4] || currNode->flag[idx][5] || currNode->flag[idx][6] || currNode->flag[idx][7])
	{
		//判断不是叶子节点
		flag = true;
	}
	else
	{
		//是叶子节点，那么保留isTarget的叶子节点的矩形范围，然后退出
		flag = false;
		//判断对应frame是否在这个叶子节点上有内容
		if (!currNode->checkNodeDataNull(isTarget))
		{
			minVList.push_back(currNode->getMinPos());
			maxVList.push_back(currNode->getMaxPos());

			//保存叶子节点的序列
			nodeList.push_back(currNode);
			currNode->setLeafFlag(leafIncr);
			leafIncr++;
		}
	}
	return flag;
}
