#include <no_comment/DTK_DealIIEntityExtraData.h>
#include <no_comment/DTK_DealIIEntityLocalMap.h>

namespace DataTransferKit {

// Partial specialization of a function is forbidden, so we delegue the work to
// some classes that can be specialized
namespace internal
{
  template <int structdim,int dim,int spacedim> struct local_map;
  template <int dim,int spacedim> struct local_map<0,dim,spacedim>;

  template <int structdim,int dim,int spacedim>
  struct local_map
  {
    static bool check_point_inclusion(
        const Teuchos::ArrayView<const double>& reference_point) 
    {
      dealii::Point<structdim> referencePoint;
      for (int d = 0; d < spacedim; ++d)
        referencePoint[d] = reference_point[d];
      // NOTE: use set parameters...
      double const epsilon = 1.0e-10;

      return dealii::GeometryInfo<dim>::is_inside_unit_cell(referencePoint, epsilon);
    }
  };



  template <int dim,int spacedim>
  struct local_map<0,dim,spacedim>
  {
    static bool check_point_inclusion(
        const Teuchos::ArrayView<const double>& reference_point) 
    {
      // Not sure what this is supposed to return for a point. So always return
      // false.
      std::ignore = reference_point;
      return false;
    }
  };
}



// helper function
template <int dim,int spacedim>
dealii::TriaIterator<dealii::CellAccessor<dim,spacedim>>
getCellIterator(Entity const & entity)
{
    Teuchos::RCP<DealIIEntityExtraData<dim,dim,spacedim>> extra_data =
        Teuchos::rcp_dynamic_cast<DealIIEntityExtraData<dim,dim,spacedim>>(entity.extraData());
    return dealii::TriaIterator<dealii::CellAccessor<dim,spacedim>>(extra_data->dealii_tria_accessor);
}



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
    auto dealii_cell_iterator = getCellIterator<dim,spacedim>(entity);
    return dealii_cell_iterator->measure();
}



template <int structdim,int dim,int spacedim>
void
DealIIEntityLocalMap<structdim,dim,spacedim>::
centroid(
    const DataTransferKit::Entity& entity,
    const Teuchos::ArrayView<double>& centroid ) const
{
    auto dealii_cell_iterator = getCellIterator<dim,spacedim>(entity);
    dealii::Point<spacedim> center_point =
        dealii_cell_iterator->center(true);
    for (int d = 0; d < spacedim; ++d)
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
    auto dealii_cell_iterator = getCellIterator<dim,spacedim>(entity);
    dealii::Point<spacedim> pointInPhysicalFrame;
    for (int d = 0; d < spacedim; ++d)
        pointInPhysicalFrame[d] = physical_point[d];

    try {
        dealii::Point<dim> pointInReferenceFrame =
            dealii_mapping->transform_real_to_unit_cell(dealii_cell_iterator, pointInPhysicalFrame);
    
        for (int d = 0; d < dim; ++d)
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
    return internal::local_map<structdim,dim,spacedim>::check_point_inclusion(reference_point);
}



template <int structdim,int dim,int spacedim>
void
DealIIEntityLocalMap<structdim,dim,spacedim>::
mapToPhysicalFrame( 
    const DataTransferKit::Entity& entity,
    const Teuchos::ArrayView<const double>& reference_point,
    const Teuchos::ArrayView<double>& physical_point ) const
{
    auto dealii_cell_iterator = getCellIterator<dim,spacedim>(entity);
    dealii::Point<dim> pointInReferenceFrame;
    for (int d = 0; d < dim; ++d)
        pointInReferenceFrame[d] = reference_point[d];

    dealii::Point<spacedim> pointInPhysicalFrame =
        dealii_mapping->transform_unit_to_real_cell(dealii_cell_iterator, pointInReferenceFrame);

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



template class DealIIEntityLocalMap<0,2,2>;
template class DealIIEntityLocalMap<0,2,3>;
template class DealIIEntityLocalMap<0,3,3>;
template class DealIIEntityLocalMap<2,2,2>;
template class DealIIEntityLocalMap<2,2,3>;
template class DealIIEntityLocalMap<3,3,3>;

} // end namespace DataTransferKit
