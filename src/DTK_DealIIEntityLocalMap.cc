#include <no_comment/DTK_DealIIEntityExtraData.h>
#include <no_comment/DTK_DealIIEntityLocalMap.h>

template <int structdim,int dim,int spacedim>
DealIIEntityLocalMap<structdim,dim,spacedim>::
DealIIEntityLocalMap(std::shared_ptr<dealii::Mapping<dim,spacedim> const> mapping)
: dealii_mapping(mapping)
{}



template <int structdim,int dim,int spacedim>
void
DealIIEntityLocalMap<structdim,dim,spacedim>::
setParameters( const Teuchos::ParameterList& parameters )
{
    std::ignore = parameters;
    // not implemented
    // probably want to read tolerances here
}



template <int structdim,int dim,int spacedim>
double
DealIIEntityLocalMap<structdim,dim,spacedim>::
measure( const DataTransferKit::Entity& entity ) const
{
    Teuchos::RCP<DealIIEntityExtraData<structdim,dim,spacedim>> extra_data
        = Teuchos::rcp_dynamic_cast<DealIIEntityExtraData<structdim,dim,spacedim>>(entity.extraData());
    auto dealii_tria_accessor = extra_data->dealii_tria_accessor;
    return dealii_tria_accessor->measure();
}



template <int structdim,int dim,int spacedim>
void
DealIIEntityLocalMap<structdim,dim,spacedim>::
centroid(
    const DataTransferKit::Entity& entity,
    const Teuchos::ArrayView<double>& centroid ) const
{
    Teuchos::RCP<DealIIEntityExtraData<structdim,dim,spacedim>> extra_data
        = Teuchos::rcp_dynamic_cast<DealIIEntityExtraData<structdim,dim,spacedim>>(entity.extraData());
    auto dealii_tria_accessor = extra_data->dealii_tria_accessor;
    dealii::Point<structdim> center_point =
        dealii_tria_accessor->center(true);
    for (int d = 0; d < structdim; ++d)
        centroid[d] = center_point[d];
}



template <int structdim,int dim,int spacedim>
bool
DealIIEntityLocalMap<structdim,dim,spacedim>::
isSafeToMapToReferenceFrame(
    const DataTransferKit::Entity& entity,
    const Teuchos::ArrayView<const double>& physical_point ) const
{
    // @Stu: don't we have a function that already check for that?
    // @Bruno: do we want to introduce a tolerance here?
    Teuchos::Tuple<double,6> bounding_box;
    entity.boundingBox(bounding_box);
    for (int d = 0; d < structdim; ++d)
        if ((physical_point[d] < bounding_box[0*dim+d])
            || (bounding_box[1*dim+d] < physical_point[d]))
        return false;
    return true;
}



template <int structdim,int dim,int spacedim>
bool
DealIIEntityLocalMap<structdim,dim,spacedim>::
mapToReferenceFrame( 
    const DataTransferKit::Entity& entity,
    const Teuchos::ArrayView<const double>& physical_point,
    const Teuchos::ArrayView<double>& reference_point ) const
{
    Teuchos::RCP<DealIIEntityExtraData<structdim,dim,spacedim>> extra_data
        = Teuchos::rcp_dynamic_cast<DealIIEntityExtraData<structdim,dim,spacedim>>(entity.extraData());
    auto dealii_cell_accessor = extra_data->dealii_tria_accessor;
    dealii::TriaIterator<dealii::CellAccessor<dim,spacedim>> dealii_tria_accessor(*dealii_cell_accessor);
    dealii::Point<spacedim> pointInPhysicalFrame;
    for (int d = 0; d < spacedim; ++d)
        pointInPhysicalFrame[d] = physical_point[d];

    try {
        dealii::Point<spacedim> pointInReferenceFrame =
            dealii_mapping->transform_real_to_unit_cell(dealii_tria_accessor, pointInPhysicalFrame);
    
        for (int d = 0; d < spacedim; ++d)
            reference_point[d] = pointInReferenceFrame[d];

        return true;
    } catch ( typename dealii::Mapping<dim,spacedim>::ExcTransformationFailed ) {
        // Maybe unnecessary...
        for (int d = 0; d < spacedim; ++d)
            reference_point[d] = std::nan("");
        return false;
    } catch ( ... ) {
        throw std::runtime_error("Unexpected exception occured!");
    }

}



template <int structdim,int dim,int spacedim>
bool
DealIIEntityLocalMap<structdim,dim,spacedim>::
checkPointInclusion( 
    const DataTransferKit::Entity& entity,
    const Teuchos::ArrayView<const double>& reference_point ) const
{
    std::ignore = entity;
    dealii::Point<structdim> referencePoint;
    for (int d = 0; d < spacedim; ++d)
        referencePoint[d] = reference_point[d];
    // NOTE: use set parameters...
    double const epsilon = 1.0e-10;
    return dealii::GeometryInfo<dim>::is_inside_unit_cell(referencePoint, epsilon);
}



template <int structdim,int dim,int spacedim>
void
DealIIEntityLocalMap<structdim,dim,spacedim>::
mapToPhysicalFrame( 
    const DataTransferKit::Entity& entity,
    const Teuchos::ArrayView<const double>& reference_point,
    const Teuchos::ArrayView<double>& physical_point ) const
{
    Teuchos::RCP<DealIIEntityExtraData<structdim,dim,spacedim>> extra_data
        = Teuchos::rcp_dynamic_cast<DealIIEntityExtraData<structdim,dim,spacedim>>(entity.extraData());
    auto dealii_cell_accessor = extra_data->dealii_tria_accessor;
    dealii::TriaIterator<dealii::CellAccessor<dim,spacedim>> dealii_tria_accessor(*dealii_cell_accessor);
    dealii::Point<spacedim> pointInReferenceFrame;
    for (int d = 0; d < spacedim; ++d)
        pointInReferenceFrame[d] = reference_point[d];

    dealii::Point<spacedim> pointInPhysicalFrame =
        dealii_mapping->transform_unit_to_real_cell(dealii_tria_accessor, pointInReferenceFrame);

    for (int d = 0; d < spacedim; ++d)
        physical_point[d] = pointInPhysicalFrame[d];
}



template <int structdim,int dim,int spacedim>
void
DealIIEntityLocalMap<structdim,dim,spacedim>::
normalAtReferencePoint( 
    const DataTransferKit::Entity& entity,
    const DataTransferKit::Entity& parent_entity,
    const Teuchos::ArrayView<const double>& reference_point,
    const Teuchos::ArrayView<double>& normal ) const
{
    std::ignore = entity;
    std::ignore = parent_entity;
    std::ignore = reference_point;
    std::ignore = normal;
    throw  std::runtime_error("DealIIEntityLocalMap::normalAtReferencePoint(...) not implemented");
}

template class DealIIEntityLocalMap<3,3,3>;
