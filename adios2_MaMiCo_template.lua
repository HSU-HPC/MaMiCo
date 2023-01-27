mmCheckVersion("")
mmCreateView("GraphEntry_1","View3DGL","::view")

mmCreateModule("BoundingBoxRenderer","::Bounding Box MD")
mmCreateModule("DistantLight","::distantlight")
mmCreateModule("SphereRenderer","::Particles")
mmCreateModule("ADIOSFlexConvert","::ADIOSFlexConvert_1")
mmCreateModule("ADIOSFlexVolume","::ADIOSFlexVolume_1")
mmCreateModule("adiosDataSource","::adiosDataSource_1")
mmCreateModule("adiosDataSource","::adiosDataSource_2")
mmCreateModule("TimeMultiplier","::TimeMultiplier_1")
mmCreateModule("RaycastVolumeRenderer","::RaycastVolumeRenderer_1")
mmCreateModule("BoundingBoxRenderer","::BoundingBoxRenderer_1")
mmCreateModule("TransferFunctionGL","::Transferfunction_volume")
mmCreateModule("TransferFunctionGL","::Transferfunction_MD")

mmCreateCall("CallRender3DGL","::view::rendering","::BoundingBoxRenderer_1::rendering")
mmCreateCall("CallRender3DGL","::Bounding Box MD::chainRendering","::Particles::rendering")
mmCreateCall("MultiParticleDataCall","::Particles::getdata","::ADIOSFlexConvert_1::mpSlot")
mmCreateCall("CallGetTransferFunctionGL","::Particles::gettransferfunction","::Transferfunction_MD::gettransferfunction")
mmCreateCall("CallLight","::Particles::lights","::distantlight::deployLightSlot")
mmCreateCall("CallADIOSData","::ADIOSFlexConvert_1::adiosSlot","::adiosDataSource_2::getdata")
mmCreateCall("CallADIOSData","::ADIOSFlexVolume_1::adiosSlot","::adiosDataSource_1::getdata")
mmCreateCall("CallRender3DGL","::TimeMultiplier_1::chainRendering","::Bounding Box MD::rendering")
mmCreateCall("CallRender3DGL","::RaycastVolumeRenderer_1::chainRendering","::TimeMultiplier_1::rendering")
mmCreateCall("VolumetricDataCall","::RaycastVolumeRenderer_1::getData","::ADIOSFlexVolume_1::volumeSlot")
mmCreateCall("CallLight","::RaycastVolumeRenderer_1::lights","::distantlight::deployLightSlot")
mmCreateCall("CallGetTransferFunctionGL","::RaycastVolumeRenderer_1::getTransferFunction","::Transferfunction_volume::gettransferfunction")
mmCreateCall("CallRender3DGL","::BoundingBoxRenderer_1::chainRendering","::RaycastVolumeRenderer_1::rendering")

mmSetParamValue("::view::camstore::settings",[=[]=])
mmSetParamValue("::view::camstore::overrideSettings",[=[false]=])
mmSetParamValue("::view::camstore::autoSaveSettings",[=[false]=])
mmSetParamValue("::view::camstore::autoLoadSettings",[=[true]=])
mmSetParamValue("::view::resetViewOnBBoxChange",[=[false]=])
mmSetParamValue("::view::showLookAt",[=[false]=])
mmSetParamValue("::view::view::showViewCube",[=[false]=])
mmSetParamValue("::view::anim::play",[=[false]=])
mmSetParamValue("::view::anim::speed",[=[4.000000]=])
mmSetParamValue("::view::anim::time",[=[1.000000]=])
mmSetParamValue("::view::backCol",[=[#000020]=])
mmSetParamValue("::view::viewKey::MoveStep",[=[0.500000]=])
mmSetParamValue("::view::viewKey::RunFactor",[=[2.000000]=])
mmSetParamValue("::view::viewKey::AngleStep",[=[90.000000]=])
mmSetParamValue("::view::viewKey::FixToWorldUp",[=[true]=])
mmSetParamValue("::view::viewKey::MouseSensitivity",[=[3.000000]=])
mmSetParamValue("::view::viewKey::RotPoint",[=[Look-At]=])
mmSetParamValue("::view::view::cubeOrientation",[=[0;0;0;1]=])
mmSetParamValue("::view::view::defaultView",[=[FACE - Front]=])
mmSetParamValue("::view::view::defaultOrientation",[=[Top]=])
mmSetParamValue("::view::cam::position",[=[-62.7262573;-140.376389;2061.86133]=])
mmSetParamValue("::view::cam::orientation",[=[0;-0;-0;1]=])
mmSetParamValue("::view::cam::projectiontype",[=[Perspective]=])
mmSetParamValue("::view::cam::nearplane",[=[1636.861328]=])
mmSetParamValue("::view::cam::farplane",[=[2366.861328]=])
mmSetParamValue("::view::cam::halfaperturedegrees",[=[28.647890]=])
mmSetParamValue("::Bounding Box MD::enableBoundingBox",[=[true]=])
mmSetParamValue("::Bounding Box MD::boundingBoxColor",[=[#ffffff]=])
mmSetParamValue("::Bounding Box MD::smoothLines",[=[true]=])
mmSetParamValue("::Bounding Box MD::enableViewCube",[=[false]=])
mmSetParamValue("::Bounding Box MD::viewCubePosition",[=[top right]=])
mmSetParamValue("::Bounding Box MD::viewCubeSize",[=[100]=])
mmSetParamValue("::distantlight::Intensity",[=[1.000000]=])
mmSetParamValue("::distantlight::Color",[=[#cccccc]=])
mmSetParamValue("::distantlight::Direction",[=[-0.25;-0.5;-0.75]=])
mmSetParamValue("::distantlight::AngularDiameter",[=[0.000000]=])
mmSetParamValue("::distantlight::EyeDirection",[=[false]=])
mmSetParamValue("::Particles::shaderMode",[=[Forward]=])
mmSetParamValue("::Particles::scaling",[=[1.000000]=])
mmSetParamValue("::Particles::forceTime",[=[false]=])
mmSetParamValue("::Particles::useLocalBbox",[=[false]=])
mmSetParamValue("::Particles::flag storage::selectedColor",[=[#ff0000]=])
mmSetParamValue("::Particles::flag storage::softSelectedColor",[=[#ff8080]=])
mmSetParamValue("::Particles::splat::alphaScaling",[=[5.000000]=])
mmSetParamValue("::Particles::splat::attenuateSubpixel",[=[false]=])
mmSetParamValue("::Particles::ssbo::staticData",[=[false]=])
mmSetParamValue("::Particles::ambient occlusion::enableLighting",[=[false]=])
mmSetParamValue("::Particles::ambient occlusion::useGsProxies",[=[false]=])
mmSetParamValue("::Particles::ambient occlusion::volumeSize",[=[128]=])
mmSetParamValue("::Particles::ambient occlusion::apex",[=[50.000000]=])
mmSetParamValue("::Particles::ambient occlusion::offset",[=[0.010000]=])
mmSetParamValue("::Particles::ambient occlusion::strength",[=[1.000000]=])
mmSetParamValue("::Particles::ambient occlusion::coneLength",[=[0.800000]=])
mmSetParamValue("::Particles::ambient occlusion::highPrecisionTexture",[=[false]=])
mmSetParamValue("::Particles::outline::width",[=[2.000000]=])
mmSetParamValue("::Particles::renderMode",[=[Simple]=])
mmSetParamValue("::ADIOSFlexConvert_1::xyzw",[=[undef]=])
mmSetParamValue("::ADIOSFlexConvert_1::xyz",[=[undef]=])
mmSetParamValue("::ADIOSFlexConvert_1::i",[=[vy]=])
mmSetParamValue("::ADIOSFlexConvert_1::box",[=[global_box]=])
mmSetParamValue("::ADIOSFlexConvert_1::x",[=[rx]=])
mmSetParamValue("::ADIOSFlexConvert_1::y",[=[ry]=])
mmSetParamValue("::ADIOSFlexConvert_1::z",[=[rz]=])
mmSetParamValue("::ADIOSFlexConvert_1::id",[=[undef]=])
mmSetParamValue("::ADIOSFlexConvert_1::direction::vx",[=[undef]=])
mmSetParamValue("::ADIOSFlexConvert_1::direction::vy",[=[undef]=])
mmSetParamValue("::ADIOSFlexConvert_1::direction::vz",[=[undef]=])
mmSetParamValue("::ADIOSFlexVolume_1::vxvyvz",[=[velocity]=])
mmSetParamValue("::ADIOSFlexVolume_1::box",[=[global_box]=])
mmSetParamValue("::ADIOSFlexVolume_1::memLayout",[=[xyz]=])
mmSetParamValue("::adiosDataSource_1::filename",[=[""]=])
mmSetParamValue("::adiosDataSource_2::filename",[=[""]=])
mmSetParamValue("::TimeMultiplier_1::multiplier",[=[200.000000]=])
mmSetParamValue("::RaycastVolumeRenderer_1::mode",[=[Integration]=])
mmSetParamValue("::RaycastVolumeRenderer_1::ray step ratio",[=[1.000000]=])
mmSetParamValue("::RaycastVolumeRenderer_1::opacity threshold",[=[1.000000]=])
mmSetParamValue("::RaycastVolumeRenderer_1::isovalue",[=[0.500000]=])
mmSetParamValue("::RaycastVolumeRenderer_1::opacity",[=[1.000000]=])
mmSetParamValue("::RaycastVolumeRenderer_1::lighting::use lighting",[=[true]=])
mmSetParamValue("::RaycastVolumeRenderer_1::lighting::ka",[=[0.100000]=])
mmSetParamValue("::RaycastVolumeRenderer_1::lighting::kd",[=[0.500000]=])
mmSetParamValue("::RaycastVolumeRenderer_1::lighting::ks",[=[0.400000]=])
mmSetParamValue("::RaycastVolumeRenderer_1::lighting::shininess",[=[10.000000]=])
mmSetParamValue("::RaycastVolumeRenderer_1::lighting::ambient color",[=[#ffffff]=])
mmSetParamValue("::RaycastVolumeRenderer_1::lighting::specular color",[=[#ffffff]=])
mmSetParamValue("::RaycastVolumeRenderer_1::lighting::material color",[=[#f2ab78]=])
mmSetParamValue("::BoundingBoxRenderer_1::enableBoundingBox",[=[true]=])
mmSetParamValue("::BoundingBoxRenderer_1::boundingBoxColor",[=[#ffffff]=])
mmSetParamValue("::BoundingBoxRenderer_1::smoothLines",[=[true]=])
mmSetParamValue("::BoundingBoxRenderer_1::enableViewCube",[=[false]=])
mmSetParamValue("::BoundingBoxRenderer_1::viewCubePosition",[=[top right]=])
mmSetParamValue("::BoundingBoxRenderer_1::viewCubeSize",[=[100]=])
mmSetParamValue("::Transferfunction_volume::TransferFunction",[=[{
  "IgnoreProjectRange": true,
  "Interpolation": "GAUSS",
  "Nodes": [
    [
      1.0,
      0.9879999756813049,
      0.0,
      0.08350914716720581,
      0.14470908045768738,
      0.012000000104308128
    ],
    [
      0.33117273449897766,
      0.824999988079071,
      0.5,
      0.125715970993042,
      0.3135363459587097,
      0.006000000052154064
    ],
    [
      0.6623454689979553,
      1.0,
      1.0,
      0.2161591649055481,
      0.46638524532318115,
      0.006000000052154064
    ]
  ],
  "TextureSize": 256,
  "ValueRange": [
    4.281161818653345e-05,
    0.37234774231910706
  ]
}]=])
mmSetParamValue("::Transferfunction_MD::TransferFunction",[=[{
  "IgnoreProjectRange": false,
  "Interpolation": "LINEAR",
  "Nodes": [
    [
      0.0,
      0.9959999918937683,
      1.0,
      1.0,
      0.0,
      0.05000000074505806
    ],
    [
      0.0,
      0.8450000286102295,
      0.8410000205039978,
      1.0,
      0.28599998354911804,
      0.05000000074505806
    ],
    [
      0.5019999742507935,
      0.5059999823570251,
      0.5019999742507935,
      1.0,
      0.5065547823905945,
      0.05000000074505806
    ],
    [
      0.8479999899864197,
      0.0,
      0.0,
      1.0,
      0.7185713648796082,
      0.05000000074505806
    ],
    [
      1.0,
      0.0,
      0.0,
      1.0,
      1.0,
      0.05000000074505806
    ]
  ],
  "TextureSize": 256,
  "ValueRange": [
    -7.0,
    7.0
  ]
}]=])

mmSetGUIVisible(true)
mmSetGUIScale(1.000000)
mmSetGUIState([=[{"AnimationEditor_State":{"animation_bounds":[0,100],"animation_file":"","animations":null,"current_frame":0,"output_prefix":"","playback_fps":30,"playing":0,"write_to_graph":false},"ConfiguratorState":{"module_list_sidebar_width":250.0,"show_module_list_sidebar":false},"GUIState":{"font_file_name":"Roboto-Regular.ttf","font_size":12.0,"global_win_background_alpha":1.0,"imgui_settings":"[Window][Configurator     F11]\nPos=843,18\nSize=1717,1087\nCollapsed=0\nDockId=0x00000004,0\n\n[Window][Parameters     F10]\nPos=0,18\nSize=841,1087\nCollapsed=0\nDockId=0x00000003,0\n\n[Window][Log Console     F9]\nPos=0,1107\nSize=2560,270\nCollapsed=0\nDockId=0x00000002,0\n\n[Window][DockSpaceViewport_11111111]\nPos=0,18\nSize=2560,1359\nCollapsed=0\n\n[Window][Debug##Default]\nPos=60,60\nSize=400,400\nCollapsed=0\n\n[Window][###fbw15804]\nPos=1080,438\nSize=400,500\nCollapsed=0\n\n[Window][###fbw17305]\nPos=1080,438\nSize=400,500\nCollapsed=0\n\n[Window][###fbw22955]\nPos=1080,438\nSize=400,500\nCollapsed=0\n\n[Window][Rename Module]\nPos=1138,651\nSize=283,74\nCollapsed=0\n\n[Window][###fbw17785]\nPos=1078,438\nSize=400,500\nCollapsed=0\n\n[Window][###fbw18471]\nPos=1080,438\nSize=400,500\nCollapsed=0\n\n[Window][###fbw19129]\nPos=1080,438\nSize=400,500\nCollapsed=0\n\n[Window][###fbw13755]\nPos=1080,438\nSize=400,500\nCollapsed=0\n\n[Window][###fbw21028]\nPos=1080,438\nSize=400,500\nCollapsed=0\n\n[Window][###fbw25215]\nPos=1080,438\nSize=400,500\nCollapsed=0\n\n[Window][###fbw27363]\nPos=1080,438\nSize=400,500\nCollapsed=0\n\n[Docking][Data]\nDockSpace     ID=0x8B93E3BD Window=0xA787BDB4 Pos=0,18 Size=2560,1359 Split=Y\n  DockNode    ID=0x00000001 Parent=0x8B93E3BD SizeRef=1920,808 Split=X\n    DockNode  ID=0x00000003 Parent=0x00000001 SizeRef=841,808 Selected=0x20D298EF\n    DockNode  ID=0x00000004 Parent=0x00000001 SizeRef=1717,808 CentralNode=1 Selected=0xB6EA6B89\n  DockNode    ID=0x00000002 Parent=0x8B93E3BD SizeRef=1920,270 Selected=0xF54B1F54\n\n","menu_visible":true,"style":2},"GraphStates":{"Project":{"Modules":{"::ADIOSFlexConvert_1":{"graph_position":[731.666748046875,390.3333435058594]},"::ADIOSFlexVolume_1":{"graph_position":[110.09874725341797,696.977294921875]},"::Bounding Box MD":{"graph_position":[-17.0,370.0]},"::BoundingBoxRenderer_1":{"graph_position":[-762.93017578125,463.3769836425781]},"::Particles":{"graph_position":[241.0,339.0]},"::RaycastVolumeRenderer_1":{"graph_position":[-509.7996826171875,438.26812744140625]},"::TimeMultiplier_1":{"graph_position":[-202.77780151367188,372.55548095703125]},"::Transferfunction_MD":{"graph_position":[857.1991577148438,531.04443359375]},"::Transferfunction_volume":{"graph_position":[199.66993713378906,911.4542846679688]},"::adiosDataSource_1":{"graph_position":[618.3334350585938,869.2223510742188]},"::adiosDataSource_2":{"graph_position":[1012.7777709960938,367.0000915527344]},"::distantlight":{"graph_position":[776.4443359375,665.33349609375]},"::view":{"graph_position":[-984.1066284179688,500.62591552734375]}},"call_coloring_map":0,"call_coloring_mode":0,"canvas_scrolling":[601.4575805664063,-114.10538482666016],"canvas_zooming":0.7938809394836426,"param_extended_mode":false,"parameter_sidebar_width":300.0,"profiling_bar_height":300.0,"project_name":"Project_1","show_call_label":true,"show_call_slots_label":false,"show_grid":false,"show_module_label":true,"show_parameter_sidebar":false,"show_profiling_bar":false,"show_slot_label":false}},"ParameterStates":{"::ADIOSFlexConvert_1::box":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::ADIOSFlexConvert_1::direction::vx":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::ADIOSFlexConvert_1::direction::vy":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::ADIOSFlexConvert_1::direction::vz":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::ADIOSFlexConvert_1::i":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::ADIOSFlexConvert_1::id":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::ADIOSFlexConvert_1::x":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::ADIOSFlexConvert_1::xyz":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::ADIOSFlexConvert_1::xyzw":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::ADIOSFlexConvert_1::y":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::ADIOSFlexConvert_1::z":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::ADIOSFlexVolume_1::box":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::ADIOSFlexVolume_1::memLayout":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::ADIOSFlexVolume_1::vxvyvz":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::Bounding Box MD::boundingBoxColor":{"gui_presentation_mode":8,"gui_read-only":false,"gui_visibility":true},"::Bounding Box MD::enableBoundingBox":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::Bounding Box MD::enableViewCube":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::Bounding Box MD::smoothLines":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::Bounding Box MD::viewCubePosition":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::Bounding Box MD::viewCubeSize":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::BoundingBoxRenderer_1::boundingBoxColor":{"gui_presentation_mode":8,"gui_read-only":false,"gui_visibility":true},"::BoundingBoxRenderer_1::enableBoundingBox":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::BoundingBoxRenderer_1::enableViewCube":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::BoundingBoxRenderer_1::smoothLines":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::BoundingBoxRenderer_1::viewCubePosition":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::BoundingBoxRenderer_1::viewCubeSize":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::Particles::ambient occlusion::apex":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":false},"::Particles::ambient occlusion::coneLength":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":false},"::Particles::ambient occlusion::enableLighting":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":false},"::Particles::ambient occlusion::highPrecisionTexture":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":false},"::Particles::ambient occlusion::offset":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":false},"::Particles::ambient occlusion::strength":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":false},"::Particles::ambient occlusion::useGsProxies":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":false},"::Particles::ambient occlusion::volumeSize":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":false},"::Particles::flag storage::selectedColor":{"gui_presentation_mode":8,"gui_read-only":false,"gui_visibility":false},"::Particles::flag storage::softSelectedColor":{"gui_presentation_mode":8,"gui_read-only":false,"gui_visibility":false},"::Particles::forceTime":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::Particles::outline::width":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":false},"::Particles::renderMode":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::Particles::scaling":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::Particles::shaderMode":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::Particles::splat::alphaScaling":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":false},"::Particles::splat::attenuateSubpixel":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":false},"::Particles::ssbo::staticData":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":false},"::Particles::useLocalBbox":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::RaycastVolumeRenderer_1::isovalue":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":false},"::RaycastVolumeRenderer_1::lighting::ambient color":{"gui_presentation_mode":8,"gui_read-only":false,"gui_visibility":true},"::RaycastVolumeRenderer_1::lighting::ka":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::RaycastVolumeRenderer_1::lighting::kd":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::RaycastVolumeRenderer_1::lighting::ks":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::RaycastVolumeRenderer_1::lighting::material color":{"gui_presentation_mode":8,"gui_read-only":false,"gui_visibility":true},"::RaycastVolumeRenderer_1::lighting::shininess":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::RaycastVolumeRenderer_1::lighting::specular color":{"gui_presentation_mode":8,"gui_read-only":false,"gui_visibility":true},"::RaycastVolumeRenderer_1::lighting::use lighting":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::RaycastVolumeRenderer_1::mode":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::RaycastVolumeRenderer_1::opacity":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":false},"::RaycastVolumeRenderer_1::opacity threshold":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::RaycastVolumeRenderer_1::ray step ratio":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::TimeMultiplier_1::multiplier":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::Transferfunction_MD::TransferFunction":{"gui_presentation_mode":32,"gui_read-only":false,"gui_visibility":true},"::Transferfunction_volume::TransferFunction":{"gui_presentation_mode":32,"gui_read-only":false,"gui_visibility":true},"::adiosDataSource_1::filename":{"gui_presentation_mode":16,"gui_read-only":false,"gui_visibility":true},"::adiosDataSource_2::filename":{"gui_presentation_mode":16,"gui_read-only":false,"gui_visibility":true},"::distantlight::AngularDiameter":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::distantlight::Color":{"gui_presentation_mode":8,"gui_read-only":false,"gui_visibility":true},"::distantlight::Direction":{"gui_presentation_mode":512,"gui_read-only":false,"gui_visibility":true},"::distantlight::EyeDirection":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::distantlight::Intensity":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::ParameterGroup::anim":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::ParameterGroup::view":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::anim::SpeedDown":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::anim::SpeedUp":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::anim::play":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::anim::speed":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::anim::time":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::anim::togglePlay":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::backCol":{"gui_presentation_mode":8,"gui_read-only":false,"gui_visibility":true},"::view::cam::farplane":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::cam::halfaperturedegrees":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::cam::nearplane":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::cam::orientation":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::cam::position":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::cam::projectiontype":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::camstore::autoLoadSettings":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::camstore::autoSaveSettings":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::camstore::overrideSettings":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::camstore::restorecam":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::camstore::settings":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::camstore::storecam":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::resetViewOnBBoxChange":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::showLookAt":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::view::cubeOrientation":{"gui_presentation_mode":2,"gui_read-only":true,"gui_visibility":false},"::view::view::defaultOrientation":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::view::defaultView":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::view::resetView":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::view::showViewCube":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::viewKey::AngleStep":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::viewKey::FixToWorldUp":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::viewKey::MouseSensitivity":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::viewKey::MoveStep":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::viewKey::RotPoint":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true},"::view::viewKey::RunFactor":{"gui_presentation_mode":2,"gui_read-only":false,"gui_visibility":true}},"WindowConfigurations":{"Animation Editor":{"win_callback":8,"win_collapsed":false,"win_flags":263168,"win_hotkey":[294,0],"win_position":[0.0,0.0],"win_reset_position":[0.0,0.0],"win_reset_size":[1600.0,800.0],"win_show":false,"win_size":[1600.0,800.0]},"Configurator":{"win_callback":6,"win_collapsed":false,"win_flags":263176,"win_hotkey":[300,0],"win_position":[843.0,18.0],"win_reset_position":[0.0,0.0],"win_reset_size":[750.0,500.0],"win_show":false,"win_size":[1717.0,1087.0]},"Hotkey Editor":{"hotkey_list":[[{"key":293,"mods":4,"name":"_hotkey_gui_exit","parent":"","parent_type":2}],[{"key":83,"mods":2,"name":"_hotkey_gui_save_project","parent":"","parent_type":2}],[{"key":76,"mods":2,"name":"_hotkey_gui_load_project","parent":"","parent_type":2}],[{"key":301,"mods":0,"name":"_hotkey_gui_menu","parent":"","parent_type":2}],[{"key":292,"mods":0,"name":"_hotkey_gui_toggle_graph_entry","parent":"","parent_type":2}],[{"key":291,"mods":0,"name":"_hotkey_gui_trigger_screenshot","parent":"","parent_type":2}],[{"key":71,"mods":2,"name":"_hotkey_gui_show-hide","parent":"","parent_type":2}],[{"key":300,"mods":0,"name":"_hotkey_gui_window_Configurator","parent":"17030246060415721360","parent_type":3}],[{"key":77,"mods":3,"name":"_hotkey_gui_configurator_module_search","parent":"","parent_type":4}],[{"key":80,"mods":3,"name":"_hotkey_gui_configurator_param_search","parent":"","parent_type":4}],[{"key":261,"mods":0,"name":"_hotkey_gui_configurator_delete_graph_entry","parent":"","parent_type":4}],[{"key":83,"mods":3,"name":"_hotkey_gui_configurator_save_project","parent":"","parent_type":4}],[{"key":82,"mods":3,"name":"_hotkey_gui_configurator_layout_graph","parent":"","parent_type":4}],[{"key":299,"mods":0,"name":"_hotkey_gui_window_Parameters","parent":"2021973547876601465","parent_type":3}],[{"key":80,"mods":2,"name":"_hotkey_gui_parameterlist_param_search","parent":"","parent_type":4}],[{"key":298,"mods":0,"name":"_hotkey_gui_window_Log Console","parent":"3923344468212134616","parent_type":3}],[{"key":297,"mods":0,"name":"_hotkey_gui_window_Transfer Function Editor","parent":"6166565350214160899","parent_type":3}],[{"key":296,"mods":0,"name":"_hotkey_gui_window_Performance Metrics","parent":"4920642770320252208","parent_type":3}],[{"key":295,"mods":0,"name":"_hotkey_gui_window_Hotkey Editor","parent":"529232018636216348","parent_type":3}],[{"key":294,"mods":0,"name":"_hotkey_gui_window_Animation Editor","parent":"12616631246881159418","parent_type":3}],[{"key":67,"mods":5,"name":"::view_camstore::storecam","parent":"::::view::camstore::storecam","parent_type":1}],[{"key":67,"mods":4,"name":"::view_camstore::restorecam","parent":"::::view::camstore::restorecam","parent_type":1}],[{"key":268,"mods":0,"name":"::view_view::resetView","parent":"::::view::view::resetView","parent_type":1}],[{"key":32,"mods":0,"name":"::view_anim::togglePlay","parent":"::::view::anim::togglePlay","parent_type":1}],[{"key":77,"mods":0,"name":"::view_anim::SpeedUp","parent":"::::view::anim::SpeedUp","parent_type":1}],[{"key":78,"mods":0,"name":"::view_anim::SpeedDown","parent":"::::view::anim::SpeedDown","parent_type":1}]],"win_callback":4,"win_collapsed":false,"win_flags":262144,"win_hotkey":[295,0],"win_position":[0.0,0.0],"win_reset_position":[0.0,0.0],"win_reset_size":[0.0,0.0],"win_show":false,"win_size":[0.0,0.0]},"Log Console":{"log_force_open":true,"log_level":2,"win_callback":7,"win_collapsed":false,"win_flags":265216,"win_hotkey":[298,0],"win_position":[0.0,1107.0],"win_reset_position":[0.0,0.0],"win_reset_size":[500.0,50.0],"win_show":true,"win_size":[2560.0,270.0]},"Parameters":{"param_extended_mode":false,"param_modules_list":[],"win_callback":1,"win_collapsed":false,"win_flags":262152,"win_hotkey":[299,0],"win_position":[0.0,18.0],"win_reset_position":[0.0,0.0],"win_reset_size":[400.0,500.0],"win_show":true,"win_size":[841.0,1087.0]},"Performance Metrics":{"fpsms_max_value_count":20,"fpsms_mode":0,"fpsms_refresh_rate":2.0,"fpsms_show_options":false,"win_callback":3,"win_collapsed":false,"win_flags":2359361,"win_hotkey":[296,0],"win_position":[0.0,0.0],"win_reset_position":[0.0,0.0],"win_reset_size":[0.0,0.0],"win_show":false,"win_size":[0.0,0.0]},"Transfer Function Editor":{"tfe_active_param":"","tfe_view_minimized":false,"tfe_view_vertical":false,"win_callback":5,"win_collapsed":false,"win_flags":262208,"win_hotkey":[297,0],"win_position":[0.0,0.0],"win_reset_position":[0.0,0.0],"win_reset_size":[0.0,0.0],"win_show":false,"win_size":[0.0,0.0]}}}]=])
