from .core import PackedPose as PackedPose, dict_to_packed as dict_to_packed, dict_to_pose as dict_to_pose, pack_result as pack_result, pose_result as pose_result, register_builtin_container_traversal as register_builtin_container_traversal, to_base64 as to_base64, to_dict as to_dict, to_packed as to_packed, to_pickle as to_pickle, to_pose as to_pose
from .pandas import register_pandas_container_traversal as register_pandas_container_traversal

def register_container_traversal(generic_func, dict_func) -> None: ...
